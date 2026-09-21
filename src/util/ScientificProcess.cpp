#include "util/ScientificProcess.h"
#include "util/TaskControl.h"
#include <stdexcept>
#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#else
#include <cerrno>
#include <csignal>
#include <fcntl.h>
#include <sys/wait.h>
#include <unistd.h>
#include <chrono>
#include <thread>
#endif

namespace
{
#ifdef _WIN32
std::wstring wide(const std::string& text)
{
    const int count = MultiByteToWideChar(CP_UTF8, MB_ERR_INVALID_CHARS, text.c_str(), -1, nullptr, 0);
    if (!count) throw std::runtime_error("Invalid UTF-8 process argument");
    std::wstring result(static_cast<std::size_t>(count), L'\0');
    MultiByteToWideChar(CP_UTF8, MB_ERR_INVALID_CHARS, text.c_str(), -1, result.data(), count);
    result.pop_back();
    return result;
}

std::wstring quoted(const std::string& text)
{
    std::wstring result = L"\"";
    std::size_t slashes = 0;
    for (wchar_t c : wide(text)) {
        if (c == L'\\') { ++slashes; continue; }
        result.append(c == L'"' ? 2 * slashes + 1 : slashes, L'\\');
        result += c;
        slashes = 0;
    }
    result.append(2 * slashes, L'\\');
    return result + L'"';
}

struct Handle
{
    HANDLE value = nullptr;
    ~Handle() { if (value && value != INVALID_HANDLE_VALUE) CloseHandle(value); }
};
#endif
}

int runScientificProcess(const std::vector<std::string>& arguments, const std::filesystem::path& logPath)
{
    if (arguments.empty() || arguments.front().empty()) throw std::runtime_error("Choose a Python interpreter");
    for (const auto& argument : arguments)
        if (argument.find('\0') != std::string::npos) throw std::runtime_error("Invalid process argument");
    atomforge::taskCheckpoint();
#ifdef _WIN32
    SECURITY_ATTRIBUTES security{sizeof(SECURITY_ATTRIBUTES), nullptr, TRUE};
    Handle log{CreateFileW(logPath.c_str(), GENERIC_WRITE, FILE_SHARE_READ, &security, CREATE_ALWAYS, FILE_ATTRIBUTE_NORMAL, nullptr)};
    if (log.value == INVALID_HANDLE_VALUE) throw std::runtime_error("Cannot create scientific process log");
    Handle input{CreateFileW(L"NUL", GENERIC_READ, FILE_SHARE_READ | FILE_SHARE_WRITE, &security, OPEN_EXISTING, 0, nullptr)};
    Handle job{CreateJobObjectW(nullptr, nullptr)};
    JOBOBJECT_EXTENDED_LIMIT_INFORMATION limits{};
    limits.BasicLimitInformation.LimitFlags = JOB_OBJECT_LIMIT_KILL_ON_JOB_CLOSE;
    if (!job.value || !SetInformationJobObject(job.value, JobObjectExtendedLimitInformation, &limits, sizeof(limits)))
        throw std::runtime_error("Cannot configure scientific process cancellation");
    std::wstring command;
    for (const auto& argument : arguments) { if (!command.empty()) command += L' '; command += quoted(argument); }
    STARTUPINFOW startup{};
    startup.cb = sizeof(startup);
    startup.dwFlags = STARTF_USESTDHANDLES;
    startup.hStdOutput = startup.hStdError = log.value;
    startup.hStdInput = input.value;
    PROCESS_INFORMATION process{};
    if (!CreateProcessW(nullptr, command.data(), nullptr, nullptr, TRUE, CREATE_NO_WINDOW | CREATE_SUSPENDED,
                        nullptr, nullptr, &startup, &process))
        throw std::runtime_error("Cannot start Python. Select an interpreter with atomforge science dependencies installed.");
    Handle processHandle{process.hProcess}, threadHandle{process.hThread};
    if (!AssignProcessToJobObject(job.value, process.hProcess)) {
        TerminateProcess(process.hProcess, 1);
        throw std::runtime_error("Cannot attach scientific process cancellation");
    }
    ResumeThread(process.hThread);
    try { while (true) {
        const DWORD status = WaitForSingleObject(process.hProcess, 100);
        if (status == WAIT_OBJECT_0) break;
        if (status == WAIT_FAILED) throw std::runtime_error("Cannot wait for scientific process");
        atomforge::taskCheckpoint();
    } } catch (...) {
        TerminateJobObject(job.value, 1);
        WaitForSingleObject(process.hProcess, INFINITE);
        throw;
    }
    DWORD code = 1;
    if (!GetExitCodeProcess(process.hProcess, &code)) throw std::runtime_error("Cannot read scientific process exit code");
    return static_cast<int>(code);
#else
    std::vector<char*> argv;
    for (const auto& argument : arguments) argv.push_back(const_cast<char*>(argument.c_str()));
    argv.push_back(nullptr);
    const int output = ::open(logPath.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0600);
    if (output < 0) throw std::runtime_error("Cannot create scientific process log");
    const pid_t child = fork();
    if (child == 0) {
        setpgid(0, 0);
        dup2(output, STDOUT_FILENO); dup2(output, STDERR_FILENO); close(output);
        const int input = ::open("/dev/null", O_RDONLY);
        if (input >= 0) { dup2(input, STDIN_FILENO); close(input); }
        execvp(argv[0], argv.data());
        _exit(127);
    }
    close(output);
    if (child < 0) throw std::runtime_error("Cannot start scientific process");
    setpgid(child, child);
    int status = 0;
    try {
        while (true) {
            const pid_t done = waitpid(child, &status, WNOHANG);
            if (done == child) break;
            if (done < 0 && errno != EINTR) throw std::runtime_error("Cannot wait for scientific process");
            atomforge::taskCheckpoint();
            std::this_thread::sleep_for(std::chrono::milliseconds(100));
        }
    } catch (...) {
        kill(-child, SIGKILL);
        while (waitpid(child, &status, 0) < 0 && errno == EINTR) {}
        throw;
    }
    return WIFEXITED(status) ? WEXITSTATUS(status) : 128;
#endif
}
