"""
Inline notebook renderer — self-contained HTML/WebGL, zero external dependencies.

Sphere impostor technique: each atom is a billboard quad; the fragment shader
ray-casts against a sphere to produce correct depth, silhouette, and Phong shading.
"""

from __future__ import annotations

import json
import math
import random
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ._structure import Structure

# Van der Waals radii in Å (fallback = 1.5)
_VDW: dict[str, float] = {
    "H": 1.20, "He": 1.40, "Li": 1.82, "Be": 1.53, "B": 1.92,
    "C": 1.70, "N": 1.55, "O": 1.52, "F": 1.47, "Ne": 1.54,
    "Na": 2.27, "Mg": 1.73, "Al": 1.84, "Si": 2.10, "P": 1.80,
    "S": 1.80, "Cl": 1.75, "Ar": 1.88, "K": 2.75, "Ca": 2.31,
    "Sc": 2.11, "Ti": 2.00, "V": 2.00, "Cr": 2.00, "Mn": 2.00,
    "Fe": 2.00, "Co": 2.00, "Ni": 1.63, "Cu": 1.40, "Zn": 1.39,
    "Ga": 1.87, "Ge": 2.11, "As": 1.85, "Se": 1.90, "Br": 1.85,
    "Kr": 2.02, "Rb": 3.03, "Sr": 2.49, "Y": 2.00, "Zr": 2.00,
    "Nb": 2.00, "Mo": 2.00, "Pd": 1.63, "Ag": 1.72, "Cd": 1.58,
    "In": 1.93, "Sn": 2.17, "Sb": 2.06, "Te": 2.06, "I": 1.98,
    "Xe": 2.16, "Cs": 3.43, "Ba": 2.68, "Pt": 1.75, "Au": 1.66,
    "Hg": 1.55, "Tl": 1.96, "Pb": 2.02, "Bi": 2.07,
    "W": 2.00, "Ta": 2.00, "Hf": 2.00, "Re": 2.00,
}

_MAX_ATOMS_RENDER = 50_000   # downsample above this for performance


def structure_to_html(s: "Structure", width: int = 620, height: int = 460) -> str:
    """Return a self-contained HTML string that renders *s* as an interactive 3D view."""

    uid = f"af{random.randint(10_000_000, 99_999_999)}"
    atoms = s.atoms
    n_total = len(atoms)

    # Downsample for performance
    step = max(1, n_total // _MAX_ATOMS_RENDER)
    atoms = atoms[::step]
    downsampled = step > 1

    # Centre the structure
    if atoms:
        cx = sum(a.x for a in atoms) / len(atoms)
        cy = sum(a.y for a in atoms) / len(atoms)
        cz = sum(a.z for a in atoms) / len(atoms)
    else:
        cx = cy = cz = 0.0

    atom_data = [
        [round(a.x - cx, 5), round(a.y - cy, 5), round(a.z - cz, 5),
         round(a.r, 4), round(a.g, 4), round(a.b, 4),
         round(_VDW.get(a.symbol, 1.5), 3)]
        for a in atoms
    ]

    # Cell box edges (12 edges from 8 corners)
    cell_lines: list = []
    if s.cell:
        av, bv, cv = s.cell
        corners = []
        for i in range(2):
            for j in range(2):
                for k in range(2):
                    corners.append([
                        round(i*av[0]+j*bv[0]+k*cv[0] - cx, 5),
                        round(i*av[1]+j*bv[1]+k*cv[1] - cy, 5),
                        round(i*av[2]+j*bv[2]+k*cv[2] - cz, 5),
                    ])
        edges = [
            (0,1),(2,3),(4,5),(6,7),  # along c
            (0,2),(1,3),(4,6),(5,7),  # along b
            (0,4),(1,5),(2,6),(3,7),  # along a
        ]
        cell_lines = [[corners[e[0]], corners[e[1]]] for e in edges]

    # Auto-scale: fit the structure to the canvas
    if atoms:
        max_ext = max(
            max(abs(a[0]) for a in atom_data),
            max(abs(a[1]) for a in atom_data),
            max(abs(a[2]) for a in atom_data),
            0.1
        )
    else:
        max_ext = 1.0

    formula_parts: dict[str, int] = {}
    for a in s.atoms:
        formula_parts[a.symbol] = formula_parts.get(a.symbol, 0) + 1
    formula = " ".join(f"{sym}{n}" for sym, n in sorted(formula_parts.items()))
    note = f" (showing 1/{step})" if downsampled else ""
    title = f"{formula}  —  {n_total} atoms{note}"

    atoms_json   = json.dumps(atom_data)
    cell_json    = json.dumps(cell_lines)
    init_scale   = round(min(width, height) * 0.38 / max(max_ext, 1.0), 3)

    # ── WebGL shaders ──────────────────────────────────────────────────────────
    vert_src = r"""
attribute vec3 aPos;       // quad corner [-1..1, -1..1, 0]
attribute vec3 aCenter;    // atom world position
attribute float aRadius;
attribute vec3 aColor;

uniform mat4 uMVP;
uniform mat4 uMV;
uniform float uScale;

varying vec3 vColor;
varying vec3 vCenterView;
varying float vRadius;
varying vec2 vQuadUV;

void main() {
    vec4 cv = uMV * vec4(aCenter, 1.0);
    vCenterView = cv.xyz;
    vRadius = aRadius * uScale;
    vColor = aColor;
    vQuadUV = aPos.xy;
    // billboard: offset in view space
    vec4 pos = vec4(cv.xyz + vec3(aPos.xy * vRadius * 1.0, 0.0), 1.0);
    gl_Position = uMVP * vec4(aCenter, 1.0);
    // override xy to billboard in clip space
    vec4 center_clip = uMVP * vec4(aCenter, 1.0);
    float aspect = 1.0;
    gl_Position = center_clip + vec4(aPos.xy * vRadius / 200.0 * center_clip.w, 0.0, 0.0);
}
"""

    frag_src = r"""
precision highp float;
varying vec3 vColor;
varying vec3 vCenterView;
varying float vRadius;
varying vec2 vQuadUV;

uniform mat4 uProj;

void main() {
    // Ray-sphere intersection in view space
    float r = length(vQuadUV);
    if (r > 1.0) discard;
    float z = sqrt(1.0 - r * r);
    vec3 normal = vec3(vQuadUV, z);

    // Phong lighting (light from upper-left in view space)
    vec3 lightDir = normalize(vec3(1.0, 1.2, 2.0));
    float diff = max(dot(normal, lightDir), 0.0);
    vec3 refl = reflect(-lightDir, normal);
    float spec = pow(max(refl.z, 0.0), 32.0);

    vec3 col = vColor * (0.25 + 0.72 * diff) + vec3(0.5) * spec * 0.35;
    gl_FragColor = vec4(col, 1.0);
}
"""

    # ── HTML template ──────────────────────────────────────────────────────────
    html = f"""
<div id="container_{uid}" style="background:#1a1a2e;border-radius:8px;
     padding:6px;display:inline-block;font-family:monospace;user-select:none;">
  <div style="color:#aaa;font-size:11px;padding:2px 6px 4px;">{title}</div>
  <canvas id="c_{uid}" width="{width}" height="{height}"
          style="display:block;cursor:grab;border-radius:4px;"></canvas>
</div>
<script>
(function(){{
  var atoms   = {atoms_json};
  var cellLines = {cell_json};
  var canvas  = document.getElementById('c_{uid}');
  var gl      = canvas.getContext('webgl') || canvas.getContext('experimental-webgl');

  // ── Fallback to Canvas 2D if WebGL unavailable ─────────────────────────────
  if (!gl) {{
    var ctx = canvas.getContext('2d');
    ctx.fillStyle = '#1a1a2e';
    ctx.fillRect(0,0,canvas.width,canvas.height);
    ctx.fillStyle = '#fff';
    ctx.font = '14px monospace';
    ctx.fillText('WebGL not available — using 2D fallback', 10, 30);
    // Simple 2D projection
    var W = canvas.width, H = canvas.height;
    var scale2d = {init_scale};
    var rotX = 0.4, rotY = 0.5;
    function project2d(x,y,z,rx,ry) {{
      var x1 = x*Math.cos(ry) + z*Math.sin(ry);
      var z1 = -x*Math.sin(ry) + z*Math.cos(ry);
      var y2 = y*Math.cos(rx) - z1*Math.sin(rx);
      var z2 = y*Math.sin(rx) + z1*Math.cos(rx);
      return [W/2 + x1*scale2d, H/2 - y2*scale2d, z2];
    }}
    function draw2d() {{
      ctx.fillStyle = '#1a1a2e';
      ctx.fillRect(0,0,W,H);
      var proj = atoms.map(function(a) {{
        var p = project2d(a[0],a[1],a[2],rotX,rotY);
        return {{px:p[0],py:p[1],pz:p[2],r:a[3],g:a[4],b:a[5],vr:a[6]}};
      }});
      proj.sort(function(a,b){{return a.pz-b.pz;}});
      proj.forEach(function(p) {{
        var sr = Math.max(2, p.vr * scale2d * 0.7);
        var ri=Math.round(p.r*255), gi=Math.round(p.g*255), bi=Math.round(p.b*255);
        var grad = ctx.createRadialGradient(p.px-sr*0.3,p.py-sr*0.3,sr*0.1,p.px,p.py,sr);
        grad.addColorStop(0,'rgba('+(ri+80)+','+(gi+80)+','+(bi+80)+',1)');
        grad.addColorStop(1,'rgba('+(ri>>1)+','+(gi>>1)+','+(bi>>1)+',1)');
        ctx.beginPath();
        ctx.arc(p.px,p.py,sr,0,Math.PI*2);
        ctx.fillStyle=grad;
        ctx.fill();
      }});
      // Draw cell box
      if (cellLines.length) {{
        ctx.strokeStyle='rgba(180,180,255,0.6)';
        ctx.lineWidth=1;
        cellLines.forEach(function(e){{
          var a=project2d(e[0][0],e[0][1],e[0][2],rotX,rotY);
          var b=project2d(e[1][0],e[1][1],e[1][2],rotX,rotY);
          ctx.beginPath(); ctx.moveTo(a[0],a[1]); ctx.lineTo(b[0],b[1]); ctx.stroke();
        }});
      }}
    }}
    var drag2d=false, lx=0, ly=0;
    canvas.addEventListener('mousedown',function(e){{drag2d=true;lx=e.offsetX;ly=e.offsetY;canvas.style.cursor='grabbing';}});
    canvas.addEventListener('mouseup',function(){{drag2d=false;canvas.style.cursor='grab';}});
    canvas.addEventListener('mousemove',function(e){{
      if(!drag2d)return;
      rotY+=(e.offsetX-lx)*0.01; rotX+=(e.offsetY-ly)*0.01;
      lx=e.offsetX; ly=e.offsetY; draw2d();
    }});
    canvas.addEventListener('wheel',function(e){{e.preventDefault();scale2d*=e.deltaY<0?1.1:0.9;draw2d();}},{{passive:false}});
    draw2d();
    return;
  }}

  // ── WebGL path ─────────────────────────────────────────────────────────────
  function compileShader(src, type) {{
    var s = gl.createShader(type);
    gl.shaderSource(s, src);
    gl.compileShader(s);
    if (!gl.getShaderParameter(s, gl.COMPILE_STATUS))
      console.error('Shader error:', gl.getShaderInfoLog(s));
    return s;
  }}
  var prog = gl.createProgram();
  gl.attachShader(prog, compileShader(`{vert_src.strip()}`, gl.VERTEX_SHADER));
  gl.attachShader(prog, compileShader(`{frag_src.strip()}`, gl.FRAGMENT_SHADER));
  gl.linkProgram(prog);
  gl.useProgram(prog);

  var aPos    = gl.getAttribLocation(prog,'aPos');
  var aCenter = gl.getAttribLocation(prog,'aCenter');
  var aRadius = gl.getAttribLocation(prog,'aRadius');
  var aColor  = gl.getAttribLocation(prog,'aColor');
  var uMVP    = gl.getUniformLocation(prog,'uMVP');
  var uMV     = gl.getUniformLocation(prog,'uMV');
  var uScale  = gl.getUniformLocation(prog,'uScale');

  // Quad vertices (2 triangles per atom billboard)
  var quadVerts = new Float32Array([-1,-1,0, 1,-1,0, 1,1,0, -1,-1,0, 1,1,0, -1,1,0]);
  var quadBuf = gl.createBuffer();
  gl.bindBuffer(gl.ARRAY_BUFFER, quadBuf);
  gl.bufferData(gl.ARRAY_BUFFER, quadVerts, gl.STATIC_DRAW);

  // Per-atom data interleaved: cx,cy,cz, radius, r,g,b  (7 floats)
  var n = atoms.length;
  var instData = new Float32Array(n * 7);
  for (var i=0; i<n; i++) {{
    var a = atoms[i];
    instData[i*7+0] = a[0]; instData[i*7+1] = a[1]; instData[i*7+2] = a[2];
    instData[i*7+3] = a[6];
    instData[i*7+4] = a[3]; instData[i*7+5] = a[4]; instData[i*7+6] = a[5];
  }}
  var instBuf = gl.createBuffer();
  gl.bindBuffer(gl.ARRAY_BUFFER, instBuf);
  gl.bufferData(gl.ARRAY_BUFFER, instData, gl.STATIC_DRAW);

  // Check for instanced arrays extension
  var ext = gl.getExtension('ANGLE_instanced_arrays');

  gl.enable(gl.DEPTH_TEST);
  gl.clearColor(0.10, 0.10, 0.18, 1.0);

  // Rotation state
  var rotX = 0.4, rotY = 0.5, scale = {init_scale};
  var drag = false, lastX = 0, lastY = 0;
  var W = canvas.width, H = canvas.height;

  function mat4Mul(a, b) {{
    var r = new Float32Array(16);
    for (var i=0;i<4;i++) for (var j=0;j<4;j++) {{
      r[i*4+j]=0; for(var k=0;k<4;k++) r[i*4+j]+=a[i*4+k]*b[k*4+j];
    }}
    return r;
  }}

  function perspective(fov, asp, near, far) {{
    var f = 1/Math.tan(fov/2);
    return new Float32Array([
      f/asp,0,0,0,
      0,f,0,0,
      0,0,(far+near)/(near-far),-1,
      0,0,(2*far*near)/(near-far),0
    ]);
  }}

  function rotMat(rx, ry) {{
    var cx=Math.cos(rx),sx=Math.sin(rx),cy=Math.cos(ry),sy=Math.sin(ry);
    return new Float32Array([
       cy,  sx*sy, -cx*sy, 0,
       0,   cx,     sx,    0,
       sy, -sx*cy,  cx*cy, 0,
       0,   0,      0,     1
    ]);
  }}

  function viewMat(dist) {{
    return new Float32Array([
      1,0,0,0,
      0,1,0,0,
      0,0,1,0,
      0,0,-dist,1
    ]);
  }}

  function draw() {{
    gl.viewport(0, 0, W, H);
    gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);

    var proj = perspective(0.7, W/H, 0.1, 5000.0);
    var view = viewMat(500);
    var rot  = rotMat(rotX, rotY);
    var mv   = mat4Mul(view, rot);
    var mvp  = mat4Mul(proj, mv);

    gl.uniformMatrix4fv(uMVP, false, mvp);
    gl.uniformMatrix4fv(uMV, false, mv);
    gl.uniform1f(uScale, scale);

    if (ext) {{
      // Instanced draw path
      gl.bindBuffer(gl.ARRAY_BUFFER, quadBuf);
      gl.enableVertexAttribArray(aPos);
      gl.vertexAttribPointer(aPos, 3, gl.FLOAT, false, 0, 0);

      gl.bindBuffer(gl.ARRAY_BUFFER, instBuf);
      var stride = 7 * 4;
      gl.enableVertexAttribArray(aCenter);
      gl.vertexAttribPointer(aCenter, 3, gl.FLOAT, false, stride, 0);
      ext.vertexAttribDivisorANGLE(aCenter, 1);

      gl.enableVertexAttribArray(aRadius);
      gl.vertexAttribPointer(aRadius, 1, gl.FLOAT, false, stride, 12);
      ext.vertexAttribDivisorANGLE(aRadius, 1);

      gl.enableVertexAttribArray(aColor);
      gl.vertexAttribPointer(aColor, 3, gl.FLOAT, false, stride, 16);
      ext.vertexAttribDivisorANGLE(aColor, 1);

      ext.drawArraysInstancedANGLE(gl.TRIANGLES, 0, 6, n);

      ext.vertexAttribDivisorANGLE(aCenter, 0);
      ext.vertexAttribDivisorANGLE(aRadius, 0);
      ext.vertexAttribDivisorANGLE(aColor, 0);
    }} else {{
      // Fallback: draw each atom individually (slower)
      for (var i = 0; i < n; i++) {{
        gl.bindBuffer(gl.ARRAY_BUFFER, quadBuf);
        gl.enableVertexAttribArray(aPos);
        gl.vertexAttribPointer(aPos, 3, gl.FLOAT, false, 0, 0);
        gl.vertexAttrib3f(aCenter, instData[i*7], instData[i*7+1], instData[i*7+2]);
        gl.vertexAttrib1f(aRadius, instData[i*7+3]);
        gl.vertexAttrib3f(aColor, instData[i*7+4], instData[i*7+5], instData[i*7+6]);
        gl.drawArrays(gl.TRIANGLES, 0, 6);
      }}
    }}
  }}

  // ── Canvas 2D cell box overlay ─────────────────────────────────────────────
  // (WebGL doesn't easily share a canvas with 2D; draw cell in same WebGL pass
  //  as line primitives if there are cell lines.)
  // For simplicity, draw cell lines as 2D canvas overlay using a second element.
  var overlay = null;
  if (cellLines.length) {{
    overlay = document.createElement('canvas');
    overlay.width = W; overlay.height = H;
    overlay.style.cssText = 'position:absolute;top:0;left:0;pointer-events:none;border-radius:4px;';
    canvas.parentNode.style.position = 'relative';
    canvas.parentNode.appendChild(overlay);
  }}

  function drawCellOverlay() {{
    if (!overlay) return;
    var ctx2 = overlay.getContext('2d');
    ctx2.clearRect(0,0,W,H);
    if (!cellLines.length) return;
    var proj = perspective(0.7, W/H, 0.1, 5000.0);
    var view = viewMat(500);
    var rot  = rotMat(rotX, rotY);
    var mv   = mat4Mul(view, rot);
    var mvp  = mat4Mul(proj, mv);
    function project(p) {{
      var x=p[0],y=p[1],z=p[2];
      var xc=mvp[0]*x+mvp[4]*y+mvp[8]*z+mvp[12];
      var yc=mvp[1]*x+mvp[5]*y+mvp[9]*z+mvp[13];
      var wc=mvp[3]*x+mvp[7]*y+mvp[11]*z+mvp[15];
      return [W/2*(1+xc/wc), H/2*(1-yc/wc)];
    }}
    ctx2.strokeStyle = 'rgba(160,180,255,0.7)';
    ctx2.lineWidth = 1.0;
    cellLines.forEach(function(e) {{
      var a = project(e[0]), b = project(e[1]);
      ctx2.beginPath(); ctx2.moveTo(a[0],a[1]); ctx2.lineTo(b[0],b[1]); ctx2.stroke();
    }});
  }}

  function redraw() {{ draw(); drawCellOverlay(); }}

  // ── Mouse & touch interaction ───────────────────────────────────────────────
  canvas.addEventListener('mousedown', function(e) {{
    drag=true; lastX=e.offsetX; lastY=e.offsetY;
    canvas.style.cursor='grabbing';
  }});
  window.addEventListener('mouseup', function() {{
    drag=false; canvas.style.cursor='grab';
  }});
  window.addEventListener('mousemove', function(e) {{
    if (!drag) return;
    var rect = canvas.getBoundingClientRect();
    var x = e.clientX - rect.left, y = e.clientY - rect.top;
    rotY += (x - lastX) * 0.01;
    rotX += (y - lastY) * 0.01;
    lastX = x; lastY = y;
    redraw();
  }});
  canvas.addEventListener('wheel', function(e) {{
    e.preventDefault();
    scale *= e.deltaY < 0 ? 1.12 : 0.89;
    redraw();
  }}, {{passive: false}});

  // Touch support
  var t0 = null, tDist0 = 0, tScale0 = 0;
  canvas.addEventListener('touchstart', function(e) {{
    e.preventDefault();
    if (e.touches.length === 1) {{
      drag=true; lastX=e.touches[0].clientX; lastY=e.touches[0].clientY;
    }} else if (e.touches.length === 2) {{
      drag=false;
      tDist0=Math.hypot(e.touches[0].clientX-e.touches[1].clientX,
                        e.touches[0].clientY-e.touches[1].clientY);
      tScale0=scale;
    }}
  }}, {{passive:false}});
  canvas.addEventListener('touchmove', function(e) {{
    e.preventDefault();
    if (e.touches.length===1 && drag) {{
      rotY+=(e.touches[0].clientX-lastX)*0.01;
      rotX+=(e.touches[0].clientY-lastY)*0.01;
      lastX=e.touches[0].clientX; lastY=e.touches[0].clientY;
      redraw();
    }} else if (e.touches.length===2) {{
      var d=Math.hypot(e.touches[0].clientX-e.touches[1].clientX,
                       e.touches[0].clientY-e.touches[1].clientY);
      scale=tScale0*(d/tDist0);
      redraw();
    }}
  }}, {{passive:false}});
  canvas.addEventListener('touchend',function(){{drag=false;}},{{passive:false}});

  redraw();
}})();
</script>
"""
    return html
