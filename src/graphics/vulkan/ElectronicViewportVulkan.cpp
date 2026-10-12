#ifdef ATOMFORGE_VULKAN
// The Vulkan side of ElectronicViewport: the same isosurface and ray-marched
// volume, drawn into a canvas (rows top first).
#include "graphics/ElectronicViewport.h"
#include "graphics/vulkan/VulkanCanvas.h"

#include "vkshaders/electronic_mesh.frag.h"
#include "vkshaders/electronic_mesh.vert.h"
#include "vkshaders/electronic_volume.frag.h"
#include "vkshaders/electronic_volume.vert.h"

#include <glm/gtc/matrix_transform.hpp>
#include <algorithm>
#include <cstddef>
#include <cstring>
#include <stdexcept>

namespace
{
struct MeshPush { glm::mat4 view; glm::vec4 clip, params; float shininess; int32_t palette; };
struct VolumePush { glm::mat4 transform; glm::vec4 dimensions, view, range; glm::ivec4 options; };
static_assert(sizeof(MeshPush)==104 && sizeof(VolumePush)==128, "push constants as the shaders");
constexpr VkShaderStageFlags kMeshStages=VK_SHADER_STAGE_VERTEX_BIT|VK_SHADER_STAGE_FRAGMENT_BIT;
const glm::vec3 kBackground(.96f,.97f,.985f);
}

struct ElectronicVulkan
{
    VkDescriptorSetLayout setLayout=VK_NULL_HANDLE;
    VkPipelineLayout meshLayout=VK_NULL_HANDLE, volumeLayout=VK_NULL_HANDLE;
    VkPipeline meshOpaque=VK_NULL_HANDLE, meshBlend=VK_NULL_HANDLE, volume=VK_NULL_HANDLE;
    VkSampler sampler=VK_NULL_HANDLE;
    bool manualFilter=false;  // R32F cannot be filtered linearly on this device
    vkr::Buffer vertices, indices;
    uint32_t indexCount=0;
    vkr::Image volumeImage;
    VkDescriptorSet volumeSet=VK_NULL_HANDLE;
    vkr::Canvas canvas;
    // The last view, redrawn by pixels().
    bool rendered=false, lastVolume=false, lastOpaque=false;
    int lastWidth=0, lastHeight=0;
    MeshPush lastMesh{};
    VolumePush lastVolumePush{};

    bool init(uint32_t stride);
    VkPipeline pipeline(const uint32_t* vertex,size_t vertexBytes,const uint32_t* fragment,size_t fragmentBytes,
                        VkPipelineLayout layout,uint32_t stride,bool mesh,bool depthWrite);
    void retireVolume();
    void draw(VkCommandBuffer cmd,VkExtent2D extent);
    ~ElectronicVulkan();
};

VkPipeline ElectronicVulkan::pipeline(const uint32_t* vertexCode,size_t vertexBytes,const uint32_t* fragmentCode,size_t fragmentBytes,
    VkPipelineLayout layout,uint32_t stride,bool mesh,bool depthWrite)
{
    vkr::Context& ctx=vkr::context();
    VkShaderModule vertex=ctx.shaderModule(vertexCode,vertexBytes), fragment=ctx.shaderModule(fragmentCode,fragmentBytes);
    VkPipelineShaderStageCreateInfo stages[2]{};
    stages[0].sType=stages[1].sType=VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO;
    stages[0].stage=VK_SHADER_STAGE_VERTEX_BIT; stages[0].module=vertex; stages[0].pName="main";
    stages[1].stage=VK_SHADER_STAGE_FRAGMENT_BIT; stages[1].module=fragment; stages[1].pName="main";
    const VkVertexInputBindingDescription binding{0,stride,VK_VERTEX_INPUT_RATE_VERTEX};
    const VkVertexInputAttributeDescription attributes[3]={{0,0,VK_FORMAT_R32G32B32_SFLOAT,0},
        {1,0,VK_FORMAT_R32G32B32_SFLOAT,12},{2,0,VK_FORMAT_R32_SFLOAT,24}};
    VkPipelineVertexInputStateCreateInfo input{VK_STRUCTURE_TYPE_PIPELINE_VERTEX_INPUT_STATE_CREATE_INFO};
    if(mesh)
    {
        input.vertexBindingDescriptionCount=1; input.pVertexBindingDescriptions=&binding;
        input.vertexAttributeDescriptionCount=3; input.pVertexAttributeDescriptions=attributes;
    }
    VkPipelineInputAssemblyStateCreateInfo assembly{VK_STRUCTURE_TYPE_PIPELINE_INPUT_ASSEMBLY_STATE_CREATE_INFO};
    assembly.topology=VK_PRIMITIVE_TOPOLOGY_TRIANGLE_LIST;
    VkPipelineViewportStateCreateInfo viewport{VK_STRUCTURE_TYPE_PIPELINE_VIEWPORT_STATE_CREATE_INFO};
    viewport.viewportCount=viewport.scissorCount=1;
    VkPipelineRasterizationStateCreateInfo raster{VK_STRUCTURE_TYPE_PIPELINE_RASTERIZATION_STATE_CREATE_INFO};
    raster.polygonMode=VK_POLYGON_MODE_FILL; raster.cullMode=VK_CULL_MODE_NONE;
    raster.frontFace=VK_FRONT_FACE_COUNTER_CLOCKWISE; raster.lineWidth=1;
    VkPipelineMultisampleStateCreateInfo multisample{VK_STRUCTURE_TYPE_PIPELINE_MULTISAMPLE_STATE_CREATE_INFO};
    multisample.rasterizationSamples=ctx.samples();
    VkPipelineDepthStencilStateCreateInfo depth{VK_STRUCTURE_TYPE_PIPELINE_DEPTH_STENCIL_STATE_CREATE_INFO};
    depth.depthTestEnable=mesh ? VK_TRUE : VK_FALSE;
    depth.depthWriteEnable=depthWrite ? VK_TRUE : VK_FALSE;
    depth.depthCompareOp=mesh ? VK_COMPARE_OP_LESS : VK_COMPARE_OP_ALWAYS;
    VkPipelineColorBlendAttachmentState blend{};
    blend.colorWriteMask=VK_COLOR_COMPONENT_R_BIT|VK_COLOR_COMPONENT_G_BIT|VK_COLOR_COMPONENT_B_BIT|VK_COLOR_COMPONENT_A_BIT;
    if(mesh)
    {
        blend.blendEnable=VK_TRUE;
        blend.srcColorBlendFactor=VK_BLEND_FACTOR_SRC_ALPHA; blend.dstColorBlendFactor=VK_BLEND_FACTOR_ONE_MINUS_SRC_ALPHA;
        blend.srcAlphaBlendFactor=VK_BLEND_FACTOR_ONE; blend.dstAlphaBlendFactor=VK_BLEND_FACTOR_ONE_MINUS_SRC_ALPHA;
        blend.colorBlendOp=blend.alphaBlendOp=VK_BLEND_OP_ADD;
    }
    VkPipelineColorBlendStateCreateInfo blendState{VK_STRUCTURE_TYPE_PIPELINE_COLOR_BLEND_STATE_CREATE_INFO};
    blendState.attachmentCount=1; blendState.pAttachments=&blend;
    const VkDynamicState dynamicStates[]={VK_DYNAMIC_STATE_VIEWPORT,VK_DYNAMIC_STATE_SCISSOR};
    VkPipelineDynamicStateCreateInfo dynamic{VK_STRUCTURE_TYPE_PIPELINE_DYNAMIC_STATE_CREATE_INFO};
    dynamic.dynamicStateCount=2; dynamic.pDynamicStates=dynamicStates;
    VkGraphicsPipelineCreateInfo info{VK_STRUCTURE_TYPE_GRAPHICS_PIPELINE_CREATE_INFO};
    info.stageCount=2; info.pStages=stages;
    info.pVertexInputState=&input; info.pInputAssemblyState=&assembly; info.pViewportState=&viewport;
    info.pRasterizationState=&raster; info.pMultisampleState=&multisample; info.pDepthStencilState=&depth;
    info.pColorBlendState=&blendState; info.pDynamicState=&dynamic;
    info.layout=layout; info.renderPass=ctx.mainRenderPass();
    VkPipeline result=VK_NULL_HANDLE;
    if(vertex && fragment) vkCreateGraphicsPipelines(ctx.device(),VK_NULL_HANDLE,1,&info,nullptr,&result);
    if(vertex) vkDestroyShaderModule(ctx.device(),vertex,nullptr);
    if(fragment) vkDestroyShaderModule(ctx.device(),fragment,nullptr);
    return result;
}

bool ElectronicVulkan::init(uint32_t stride)
{
    vkr::Context& ctx=vkr::context();
    if(!ctx.ready() || !vkr::sceneRenderer().ready()) return false;
    const VkDevice device=ctx.device();
    VkFormatProperties format{};
    vkGetPhysicalDeviceFormatProperties(ctx.physicalDevice(),VK_FORMAT_R32_SFLOAT,&format);
    manualFilter=!(format.optimalTilingFeatures & VK_FORMAT_FEATURE_SAMPLED_IMAGE_FILTER_LINEAR_BIT);
    VkSamplerCreateInfo samplerInfo{VK_STRUCTURE_TYPE_SAMPLER_CREATE_INFO};
    samplerInfo.magFilter=samplerInfo.minFilter=manualFilter ? VK_FILTER_NEAREST : VK_FILTER_LINEAR;
    samplerInfo.addressModeU=samplerInfo.addressModeV=samplerInfo.addressModeW=VK_SAMPLER_ADDRESS_MODE_CLAMP_TO_EDGE;
    vkCreateSampler(device,&samplerInfo,nullptr,&sampler);
    VkDescriptorSetLayoutBinding binding{0,VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER,1,VK_SHADER_STAGE_FRAGMENT_BIT,nullptr};
    VkDescriptorSetLayoutCreateInfo setInfo{VK_STRUCTURE_TYPE_DESCRIPTOR_SET_LAYOUT_CREATE_INFO};
    setInfo.bindingCount=1; setInfo.pBindings=&binding;
    vkCreateDescriptorSetLayout(device,&setInfo,nullptr,&setLayout);
    const VkPushConstantRange meshRange{kMeshStages,0,sizeof(MeshPush)}, volumeRange{VK_SHADER_STAGE_FRAGMENT_BIT,0,sizeof(VolumePush)};
    VkPipelineLayoutCreateInfo layoutInfo{VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO};
    layoutInfo.pushConstantRangeCount=1; layoutInfo.pPushConstantRanges=&meshRange;
    vkCreatePipelineLayout(device,&layoutInfo,nullptr,&meshLayout);
    layoutInfo.setLayoutCount=1; layoutInfo.pSetLayouts=&setLayout; layoutInfo.pPushConstantRanges=&volumeRange;
    vkCreatePipelineLayout(device,&layoutInfo,nullptr,&volumeLayout);
    if(!sampler || !setLayout || !meshLayout || !volumeLayout) return false;
    meshOpaque=pipeline(kSpv_electronic_mesh_vert,sizeof(kSpv_electronic_mesh_vert),kSpv_electronic_mesh_frag,sizeof(kSpv_electronic_mesh_frag),
                        meshLayout,stride,true,true);
    meshBlend=pipeline(kSpv_electronic_mesh_vert,sizeof(kSpv_electronic_mesh_vert),kSpv_electronic_mesh_frag,sizeof(kSpv_electronic_mesh_frag),
                       meshLayout,stride,true,false);
    volume=pipeline(kSpv_electronic_volume_vert,sizeof(kSpv_electronic_volume_vert),kSpv_electronic_volume_frag,sizeof(kSpv_electronic_volume_frag),
                    volumeLayout,stride,false,false);
    return meshOpaque && meshBlend && volume;
}

void ElectronicVulkan::retireVolume()
{
    vkr::Context& ctx=vkr::context();
    ctx.destroyImage(volumeImage);
    if(volumeSet)
    {
        const VkDevice device=ctx.device();
        const VkDescriptorPool pool=ctx.descriptorPool();
        const VkDescriptorSet set=volumeSet;
        ctx.destroyLater([device,pool,set] { vkFreeDescriptorSets(device,pool,1,&set); });
        volumeSet=VK_NULL_HANDLE;
    }
}

ElectronicVulkan::~ElectronicVulkan()
{
    vkr::Context& ctx=vkr::context();
    ctx.destroyBuffer(vertices); ctx.destroyBuffer(indices);
    retireVolume();
    const VkDevice device=ctx.device();
    const VkPipeline pipelines[]={meshOpaque,meshBlend,volume};
    const VkPipelineLayout layouts[]={meshLayout,volumeLayout};
    const VkDescriptorSetLayout set=setLayout;
    const VkSampler doomedSampler=sampler;
    ctx.destroyLater([device,pipelines,layouts,set,doomedSampler] {
        for(auto p:pipelines) if(p) vkDestroyPipeline(device,p,nullptr);
        for(auto l:layouts) if(l) vkDestroyPipelineLayout(device,l,nullptr);
        if(set) vkDestroyDescriptorSetLayout(device,set,nullptr);
        if(doomedSampler) vkDestroySampler(device,doomedSampler,nullptr);
    });
}

void ElectronicVulkan::draw(VkCommandBuffer cmd,VkExtent2D extent)
{
    const VkViewport viewport{0,0,static_cast<float>(extent.width),static_cast<float>(extent.height),0,1};
    const VkRect2D scissor{{0,0},extent};
    vkCmdSetViewport(cmd,0,1,&viewport); vkCmdSetScissor(cmd,0,1,&scissor);
    if(lastVolume)
    {
        if(!volumeSet) return;
        vkCmdBindPipeline(cmd,VK_PIPELINE_BIND_POINT_GRAPHICS,volume);
        vkCmdBindDescriptorSets(cmd,VK_PIPELINE_BIND_POINT_GRAPHICS,volumeLayout,0,1,&volumeSet,0,nullptr);
        vkCmdPushConstants(cmd,volumeLayout,VK_SHADER_STAGE_FRAGMENT_BIT,0,sizeof(VolumePush),&lastVolumePush);
        vkCmdDraw(cmd,3,1,0,0);
        return;
    }
    if(!vertices.buffer || !indices.buffer || !indexCount) return;
    vkCmdBindPipeline(cmd,VK_PIPELINE_BIND_POINT_GRAPHICS,lastOpaque ? meshOpaque : meshBlend);
    const VkDeviceSize offset=0;
    vkCmdBindVertexBuffers(cmd,0,1,&vertices.buffer,&offset);
    vkCmdBindIndexBuffer(cmd,indices.buffer,0,VK_INDEX_TYPE_UINT32);
    vkCmdPushConstants(cmd,meshLayout,kMeshStages,0,sizeof(MeshPush),&lastMesh);
    vkCmdDrawIndexed(cmd,indexCount,1,0,0,0);
}

namespace
{
// The renderer, created once the Vulkan context is ready (null before).
ElectronicVulkan* vulkan(std::shared_ptr<ElectronicVulkan>& vk,uint32_t stride)
{
    if(!vk)
    {
        auto created=std::make_shared<ElectronicVulkan>();
        if(!created->init(stride)) return nullptr;
        vk=created;
    }
    return vk.get();
}

ElectronicViewport::Texture present(ElectronicVulkan& vk)
{
    vkr::Context& ctx=vkr::context();
    vk.rendered=true;
    if(ctx.inFrame())
    {
        if(!vk.canvas.begin(vk.lastWidth,vk.lastHeight,kBackground)) return 0;
        vk.draw(vk.canvas.commandBuffer(),vk.canvas.extent());
        vk.canvas.end();
    }
    else ctx.immediate([&](VkCommandBuffer cmd) {
        if(!vk.canvas.begin(vk.lastWidth,vk.lastHeight,kBackground,cmd)) return;
        vk.draw(cmd,vk.canvas.extent());
        vk.canvas.end();
    });
    return static_cast<ElectronicViewport::Texture>(vk.canvas.texture());
}
}

ElectronicViewport::Texture ElectronicViewport::renderVulkan(int width,int height,float yaw,float pitch,float zoom,glm::vec2 pan,
    float opacity,float low,float high,int palette,float specular,float shininess)
{
    static_assert(sizeof(Vertex)==28 && offsetof(Vertex,normal)==12 && offsetof(Vertex,scalar)==24, "vertex layout as the pipeline");
    ElectronicVulkan* vk=vulkan(m_vk,sizeof(Vertex));
    if(!vk) return 0;
    vkr::Context& ctx=vkr::context();
    if(m_dirty)
    {
        ctx.destroyBuffer(vk->vertices);
        if(!m_vertices.empty()) vk->vertices=ctx.createDeviceBuffer(m_vertices.data(),m_vertices.size()*sizeof(Vertex),VK_BUFFER_USAGE_VERTEX_BUFFER_BIT);
    }
    if(sortTriangles(yaw,pitch))
    {
        // A new buffer: the previous one may still be read by a frame in flight.
        ctx.destroyBuffer(vk->indices); vk->indexCount=0;
        if(!m_indices.empty())
        {
            vk->indices=ctx.createBuffer(m_indices.size()*sizeof(unsigned int),VK_BUFFER_USAGE_INDEX_BUFFER_BIT,true);
            if(!vk->indices.mapped) throw std::runtime_error("Electronic viewport index buffer could not be allocated");
            std::memcpy(vk->indices.mapped,m_indices.data(),m_indices.size()*sizeof(unsigned int));
            vk->indexCount=static_cast<uint32_t>(m_indices.size());
        }
        m_dirty=false;
    }
    if(!m_vertices.empty() && !vk->vertices.buffer) throw std::runtime_error("Electronic viewport vertex buffer could not be allocated");
    const auto rotation=glm::rotate(glm::rotate(glm::mat4(1),pitch,glm::vec3(1,0,0)),yaw,glm::vec3(0,1,0));
    const float aspect=static_cast<float>(width)/height, extent=1.25f/std::max(.1f,zoom);
    const auto projection=vkr::clipCorrection(true)*glm::ortho(-extent*aspect,extent*aspect,-extent,extent,.1f,10.0f);
    MeshPush& push=vk->lastMesh;
    push.view=glm::translate(glm::mat4(1),glm::vec3(pan.x,pan.y,-3))*rotation;
    push.clip=glm::vec4(projection[0][0],projection[1][1],projection[2][2],projection[3][2]);
    push.params=glm::vec4(low,high,opacity,specular);
    push.shininess=shininess; push.palette=palette;
    vk->lastOpaque=opacity>=.999f; vk->lastVolume=false;
    vk->lastWidth=width; vk->lastHeight=height;
    m_width=width; m_height=height;
    return present(*vk);
}

void ElectronicViewport::setVolumeVulkan(const std::vector<float>& values)
{
    ElectronicVulkan* vk=vulkan(m_vk,sizeof(Vertex));
    if(!vk) return;
    vkr::Context& ctx=vkr::context();
    const VkDevice device=ctx.device();
    vkr::Image image;
    image.format=VK_FORMAT_R32_SFLOAT;
    VkImageCreateInfo info{VK_STRUCTURE_TYPE_IMAGE_CREATE_INFO};
    info.imageType=VK_IMAGE_TYPE_3D; info.format=image.format;
    info.extent={static_cast<uint32_t>(m_volumeShape.x),static_cast<uint32_t>(m_volumeShape.y),static_cast<uint32_t>(m_volumeShape.z)};
    info.mipLevels=info.arrayLayers=1; info.samples=VK_SAMPLE_COUNT_1_BIT; info.tiling=VK_IMAGE_TILING_OPTIMAL;
    info.usage=VK_IMAGE_USAGE_SAMPLED_BIT|VK_IMAGE_USAGE_TRANSFER_DST_BIT;
    info.sharingMode=VK_SHARING_MODE_EXCLUSIVE; info.initialLayout=VK_IMAGE_LAYOUT_UNDEFINED;
    if(vkCreateImage(device,&info,nullptr,&image.image)!=VK_SUCCESS) throw std::runtime_error("Volume texture could not be created");
    VkMemoryRequirements requirements;
    vkGetImageMemoryRequirements(device,image.image,&requirements);
    VkMemoryAllocateInfo allocate{VK_STRUCTURE_TYPE_MEMORY_ALLOCATE_INFO};
    allocate.allocationSize=requirements.size;
    allocate.memoryTypeIndex=ctx.memoryType(requirements.memoryTypeBits,VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT);
    if(allocate.memoryTypeIndex==UINT32_MAX || vkAllocateMemory(device,&allocate,nullptr,&image.memory)!=VK_SUCCESS)
    {
        vkDestroyImage(device,image.image,nullptr);
        throw std::runtime_error("Not enough GPU memory for the volume texture");
    }
    vkBindImageMemory(device,image.image,image.memory,0);
    VkImageViewCreateInfo viewInfo{VK_STRUCTURE_TYPE_IMAGE_VIEW_CREATE_INFO};
    viewInfo.image=image.image; viewInfo.viewType=VK_IMAGE_VIEW_TYPE_3D; viewInfo.format=image.format;
    viewInfo.subresourceRange={VK_IMAGE_ASPECT_COLOR_BIT,0,1,0,1};
    vkCreateImageView(device,&viewInfo,nullptr,&image.view);
    vkr::Buffer staging=ctx.createBuffer(values.size()*sizeof(float),VK_BUFFER_USAGE_TRANSFER_SRC_BIT,true);
    if(!image.view || !staging.mapped)
    {
        ctx.destroyBuffer(staging); ctx.destroyImage(image);
        throw std::runtime_error("Not enough memory to upload the volume texture");
    }
    std::memcpy(staging.mapped,values.data(),values.size()*sizeof(float));
    ctx.immediate([&](VkCommandBuffer cmd) {
        VkImageMemoryBarrier barrier{VK_STRUCTURE_TYPE_IMAGE_MEMORY_BARRIER};
        barrier.srcAccessMask=0; barrier.dstAccessMask=VK_ACCESS_TRANSFER_WRITE_BIT;
        barrier.oldLayout=VK_IMAGE_LAYOUT_UNDEFINED; barrier.newLayout=VK_IMAGE_LAYOUT_TRANSFER_DST_OPTIMAL;
        barrier.srcQueueFamilyIndex=barrier.dstQueueFamilyIndex=VK_QUEUE_FAMILY_IGNORED;
        barrier.image=image.image; barrier.subresourceRange={VK_IMAGE_ASPECT_COLOR_BIT,0,1,0,1};
        vkCmdPipelineBarrier(cmd,VK_PIPELINE_STAGE_TOP_OF_PIPE_BIT,VK_PIPELINE_STAGE_TRANSFER_BIT,0,0,nullptr,0,nullptr,1,&barrier);
        VkBufferImageCopy region{};
        region.imageSubresource={VK_IMAGE_ASPECT_COLOR_BIT,0,0,1};
        region.imageExtent=info.extent;
        vkCmdCopyBufferToImage(cmd,staging.buffer,image.image,VK_IMAGE_LAYOUT_TRANSFER_DST_OPTIMAL,1,&region);
        barrier.srcAccessMask=VK_ACCESS_TRANSFER_WRITE_BIT; barrier.dstAccessMask=VK_ACCESS_SHADER_READ_BIT;
        barrier.oldLayout=VK_IMAGE_LAYOUT_TRANSFER_DST_OPTIMAL; barrier.newLayout=VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL;
        vkCmdPipelineBarrier(cmd,VK_PIPELINE_STAGE_TRANSFER_BIT,VK_PIPELINE_STAGE_FRAGMENT_SHADER_BIT,0,0,nullptr,0,nullptr,1,&barrier);
    });
    ctx.destroyBuffer(staging);
    VkDescriptorSet set=VK_NULL_HANDLE;
    VkDescriptorSetAllocateInfo setInfo{VK_STRUCTURE_TYPE_DESCRIPTOR_SET_ALLOCATE_INFO};
    setInfo.descriptorPool=ctx.descriptorPool(); setInfo.descriptorSetCount=1; setInfo.pSetLayouts=&vk->setLayout;
    if(vkAllocateDescriptorSets(device,&setInfo,&set)!=VK_SUCCESS)
    {
        ctx.destroyImage(image);
        throw std::runtime_error("Volume texture descriptor could not be allocated");
    }
    const VkDescriptorImageInfo imageInfo{vk->sampler,image.view,VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL};
    VkWriteDescriptorSet write{VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET};
    write.dstSet=set; write.dstBinding=0; write.descriptorCount=1;
    write.descriptorType=VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER; write.pImageInfo=&imageInfo;
    vkUpdateDescriptorSets(device,1,&write,0,nullptr);
    // The previous texture may still be read by a frame in flight.
    vk->retireVolume();
    vk->volumeImage=image; vk->volumeSet=set;
}

ElectronicViewport::Texture ElectronicViewport::renderVolumeVulkan(int width,int height,float yaw,float pitch,float zoom,glm::vec2 pan,
    float opacity,float threshold,float low,float high,int palette)
{
    ElectronicVulkan* vk=vulkan(m_vk,sizeof(Vertex));
    if(!vk) return 0;
    const auto inverseRotation=glm::transpose(glm::mat3(glm::rotate(glm::rotate(glm::mat4(1),pitch,glm::vec3(1,0,0)),yaw,glm::vec3(0,1,0))));
    VolumePush& push=vk->lastVolumePush;
    push.transform=glm::mat4(m_inverseCell*inverseRotation);
    push.dimensions=glm::vec4(glm::vec3(m_volumeShape),static_cast<float>(width)/height);
    push.view=glm::vec4(pan.x,pan.y,1.25f/std::max(.1f,zoom),opacity);
    push.range=glm::vec4(threshold,low,high,0);
    push.options=glm::ivec4(m_raySteps,palette,vk->manualFilter ? 1 : 0,0);
    vk->lastVolume=true;
    vk->lastWidth=width; vk->lastHeight=height;
    m_width=width; m_height=height;
    return present(*vk);
}

std::vector<unsigned char> ElectronicViewport::pixelsVulkan() const
{
    if(!m_vk || !m_vk->rendered || m_width<=0 || m_height<=0) throw std::runtime_error("Render a view before exporting");
    ElectronicVulkan& vk=*m_vk;
    vkr::Context& ctx=vkr::context();
    const std::size_t bytes=static_cast<std::size_t>(vk.lastWidth)*vk.lastHeight*4;
    vkr::Buffer readback=ctx.createBuffer(bytes,VK_BUFFER_USAGE_TRANSFER_DST_BIT,true);
    if(!readback.mapped) throw std::runtime_error("Not enough memory to export the view");
    // Drawn again at once: the canvas image of this frame is not rendered yet.
    vkr::Canvas canvas;
    bool drawn=false;
    ctx.immediate([&](VkCommandBuffer cmd) {
        if(!canvas.begin(vk.lastWidth,vk.lastHeight,kBackground,cmd)) return;
        vk.draw(cmd,canvas.extent());
        canvas.end();
        canvas.copyTo(cmd,readback);
        drawn=true;
    });
    std::vector<unsigned char> result;
    if(drawn)
    {
        result.resize(bytes);
        std::memcpy(result.data(),readback.mapped,bytes);
        if(ctx.colorFormat()==VK_FORMAT_B8G8R8A8_UNORM || ctx.colorFormat()==VK_FORMAT_B8G8R8A8_SRGB)
            for(std::size_t i=0;i<bytes;i+=4) std::swap(result[i],result[i+2]);
        for(std::size_t i=3;i<bytes;i+=4) result[i]=255;
    }
    ctx.destroyBuffer(readback);
    if(!drawn) throw std::runtime_error("Failed to draw the view for export");
    return result;
}
#endif
