#!/usr/bin/env python3
import os
import sys
import argparse
import shutil
import subprocess
import platform
from colorama import Fore, Style, init

def check_dependencies():
    """Check for required dependencies and return missing ones."""
    print(f"{Fore.CYAN}Checking dependencies...{Style.RESET_ALL}")
    
    is_linux = sys.platform.startswith("linux")
    is_macos = sys.platform == "darwin"
    
    # Define dependencies based on platform
    if is_linux:
        dependencies = {
            "cmake": {"cmd": "cmake", "pkg": "cmake"},
            "g++": {"cmd": "g++", "pkg": "build-essential"},
            "clang++": {"cmd": "clang++", "pkg": "clang"},
            "qmake": {"cmd": "qmake", "pkg": "qtbase5-dev qt5-qmake qtbase5-dev-tools"},
            "ninja": {"cmd": "ninja", "pkg": "ninja-build"},
            "rdkit": {"file": "/usr/include/rdkit/RDGeneral/RDLog.h", "pkg": "librdkit-dev"},
            "eigen": {"file": "/usr/include/eigen3/Eigen/Core", "pkg": "libeigen3-dev"},
            "tbb": {"file": "/usr/include/tbb/tbb.h", "pkg": "libtbb-dev"},
            "rapidjson": {"file": "/usr/include/rapidjson/rapidjson.h", "pkg": "rapidjson-dev"},
            "cxxopts": {"file": "/usr/include/cxxopts.hpp", "pkg": "libcxxopts-dev"}
        }
    elif is_macos:
        dependencies = {
            "cmake": {"cmd": "cmake", "pkg": "cmake"},
            "clang++": {"cmd": "clang++", "pkg": "xcode-select --install"},
            "qmake": {"cmd": "qmake", "pkg": "qt@5"},
            "ninja": {"cmd": "ninja", "pkg": "ninja"},
            "rdkit": {"file": "/usr/local/include/rdkit/RDGeneral/RDLog.h", "pkg": "rdkit"},
            "eigen": {"file": "/usr/local/include/eigen3/Eigen/Core", "pkg": "eigen"},
            "tbb": {"file": "/usr/local/include/tbb/tbb.h", "pkg": "tbb"},
            "rapidjson": {"file": "/usr/local/include/rapidjson/rapidjson.h", "pkg": "rapidjson"},
            "cxxopts": {"file": "/usr/local/include/cxxopts.hpp", "pkg": "cxxopts"}
        }
    else:
        print(f"{Fore.YELLOW}Dependency checking not implemented for this platform.{Style.RESET_ALL}")
        return []
    
    missing = []
    
    for name, attrs in dependencies.items():
        found = False
        if "cmd" in attrs:
            cmd_path = shutil.which(attrs["cmd"])
            if cmd_path:
                if name == "cmake" or name == "ninja":
                    version_info = "unknown"
                    try:
                        result = subprocess.run([cmd_path, "--version"], 
                                            stdout=subprocess.PIPE, 
                                            stderr=subprocess.STDOUT, 
                                            text=True, 
                                            timeout=5)
                        version_info = result.stdout.strip().split('\n')[0]
                    except:
                        pass
                    print(f"{Fore.GREEN}✓ Found {name}: {version_info}{Style.RESET_ALL}")
                else:
                    print(f"{Fore.GREEN}✓ Found {name}{Style.RESET_ALL}")
                found = True
        
        if not found and "file" in attrs:
            if os.path.exists(attrs["file"]):
                print(f"{Fore.GREEN}✓ Found {name}{Style.RESET_ALL}")
                found = True
        
        if not found:
            print(f"{Fore.RED}✗ Missing {name}{Style.RESET_ALL}")
            missing.append(attrs["pkg"])
    
    return missing

def display_install_instructions(missing_deps):
    """Display instructions for installing missing dependencies."""
    if not missing_deps:
        return
    
    is_linux = sys.platform.startswith("linux")
    is_macos = sys.platform == "darwin"
    
    print(f"\n{Fore.RED}Missing dependencies found!{Style.RESET_ALL}")
    
    if is_linux:
        print(f"\n{Fore.YELLOW}To install missing dependencies on Ubuntu/Debian:{Style.RESET_ALL}")
        print(f"sudo apt install {' '.join(missing_deps)}")
    elif is_macos:
        print(f"\n{Fore.YELLOW}To install missing dependencies on macOS:{Style.RESET_ALL}")
        print(f"brew install {' '.join(missing_deps)}")
    
    print(f"\n{Fore.YELLOW}After installing the dependencies, run this script again.{Style.RESET_ALL}")
    
    if get_user_choice("Continue anyway?", ["Yes", "No"], "No") == "No":
        sys.exit(1)

def show_help_for_option(option_type):
    """Display help text for different option types."""
    help_texts = {
        "optimization": """
Optimization Levels:
  -O0     - No optimization, fastest compile time, best for debugging
  -O1     - Basic optimization, reasonable compile time and better performance than O0
  -O2     - More aggressive optimization, good balance between performance and compile time
  -O3     - Maximum optimization, slower compile time but best runtime performance
  -Ofast  - Beyond O3, may break strict standard compliance for maximum speed

Choose O0 for debugging, O2 for development, and O3/Ofast for production releases.
        """,
        
        "cpu": """
CPU Instruction Sets:
  -march=native    - Optimizes for your exact CPU model
  -mavx/-mavx2     - Uses AVX/AVX2 SIMD instructions (faster vector operations)
  -mavx512*        - Uses AVX-512 instructions (available on newer Intel CPUs)
  -msse4.2         - Uses SSE4.2 instructions (older SIMD standard)
  
Native is the best choice for single-machine deployment.
Specific instruction sets are better for distributable binaries.
        """,
        
        "advanced": """
Advanced Optimizations:
  -flto          - Link Time Optimization combines optimization across all files
  -fprofile-*    - Profile-Guided Optimization uses runtime data to improve optimization
  -polly         - Polyhedral loop optimization (Clang only)
  
LTO adds compile time but improves performance significantly.
PGO requires multi-step build but provides best performance.
        """,
        
        "debug": """
Debug Information:
  -g0    - No debug information
  -g1    - Minimal debug information 
  -g2    - Standard debug information
  -g3    - Maximum debug information
  
Higher numbers increase binary size but improve debuggability.
        """,
        
        "gui": """
GUI Option:
  The GUI provides a graphical interface for DescriptorFactory.
  It requires Qt libraries to be installed.
  
  Enable if you want a visual interface for managing molecular descriptors.
  Disable for server deployments or command-line only usage.
        """,
        
        "build_system": """
Build Systems:
  CMake+Ninja   - Fastest build system, good parallel builds
  CMake+Make    - Traditional build system, widely available
  Make          - Direct Makefile build without CMake
  
Ninja is recommended for development (faster builds).
CMake+Make is the most portable option.
        """,
    }
    
    if option_type in help_texts:
        print(f"\n{Fore.CYAN}=== Help: {option_type.replace('_', ' ').title()} ==={Style.RESET_ALL}")
        print(help_texts[option_type])
        input(f"{Fore.YELLOW}Press Enter to continue...{Style.RESET_ALL}")
    else:
        print(f"{Fore.RED}No help available for {option_type}{Style.RESET_ALL}")
        input(f"{Fore.YELLOW}Press Enter to continue...{Style.RESET_ALL}")

def get_user_choice(prompt, options, default=None, help_key=None):
    options_lower = [opt.lower() for opt in options]
    default_lower = default.lower() if default else None

    if default and default_lower not in options_lower:
        default = None
        default_lower = None

    prompt_options = "/".join([f"{Style.BRIGHT}{opt.upper()}{Style.RESET_ALL}" if opt.lower() == default_lower else opt for opt in options])
    full_prompt = f"{Fore.CYAN}?{Style.RESET_ALL} {prompt} ({prompt_options}) or 'q' to exit"
    if help_key:
        full_prompt += " or '?' for help"
    full_prompt += ": "

    while True:
        choice = input(full_prompt).strip().lower()
        if choice == 'q':
            print(f"{Fore.YELLOW}Exiting...{Style.RESET_ALL}")
            sys.exit(0)
        elif choice == '?' and help_key:
            show_help_for_option(help_key)
            continue
        if not choice and default:
            return default
        if choice in options_lower:
            return options[options_lower.index(choice)]
        print(f"{Fore.RED}Invalid choice.{Style.RESET_ALL} Please enter one of: {', '.join(options)}, or 'q' to exit")

def get_compiler_info():
    compilers = {}
    
    for compiler_name in ["clang++", "g++", "c++"]:
        compiler_path = shutil.which(compiler_name)
        if compiler_path:
            try:
                result = subprocess.run([compiler_path, "--version"],
                                     stdout=subprocess.PIPE,
                                     stderr=subprocess.STDOUT,
                                     text=True,
                                     timeout=5)
                
                version_info = result.stdout.strip()
                first_line = version_info.split('\n')[0]

                compiler_type = "unknown"
                if "clang" in version_info.lower():
                    compiler_type = "clang"
                elif "gcc" in version_info.lower() or "g++" in version_info.lower():
                    compiler_type = "gcc"

                if version_info and compiler_type != "unknown":
                    compilers[compiler_name] = {
                        "version": first_line,
                        "path": compiler_path,
                        "type": compiler_type
                    }
                    print(f"{Fore.GREEN}Found {compiler_type.upper()} ({compiler_name}): {first_line}{Style.RESET_ALL}")
            except Exception:
                pass
    
    if not compilers:
        print(f"{Fore.RED}No suitable C++ compiler found in PATH.{Style.RESET_ALL}")
    
    return compilers

def detect_cpu_capabilities():
    capabilities = []
    machine = platform.machine()
    
    try:
        if sys.platform == "linux":
            with open("/proc/cpuinfo", "r") as f:
                cpu_info = f.read()
            flags_match = None
            for line in cpu_info.splitlines():
                if line.startswith("flags"):
                    flags_match = line.split(":", 1)[1].strip().split()
                    break
                    
            if flags_match:
                flags_lower = [f.lower() for f in flags_match]
                for feature in ["avx", "avx2", "avx512f", "fma", "sse4_2"]:
                    if feature in flags_lower:
                        if feature == "avx512f":
                            if "avx512" not in capabilities:
                                capabilities.append("avx512")
                        else:
                            if feature not in capabilities:
                                capabilities.append(feature)
                
        elif sys.platform == "darwin":
            result = subprocess.run(["sysctl", "machdep.cpu.features"], 
                                 stdout=subprocess.PIPE, text=True, check=False)
            if result.returncode == 0:
                features = result.stdout.split(":", 1)[1].strip().upper().split()
                for feature in ["AVX", "AVX2", "AVX512F", "FMA", "SSE4_2"]:
                    if feature in features:
                        cap_name = feature.lower()
                        if cap_name == "avx512f":
                            if "avx512" not in capabilities:
                                capabilities.append("avx512")
                        else:
                            if cap_name not in capabilities:
                                capabilities.append(cap_name)
                                
    except Exception:
        pass
        
    return capabilities

def get_optimization_flags(compiler_type, capabilities):
    flag_options = {}
    
    # Basic optimization levels
    flag_options["Optimization Level"] = {
        "O0": "No optimization",
        "O1": "Basic optimization", 
        "O2": "Recommended optimization",
        "O3": "Aggressive optimization",
        "Ofast": "Maximum optimization (may affect precision)"
    }
    
    # CPU-specific instructions
    cpu_flags = {}
    if "avx512" in capabilities:
        cpu_flags["avx512"] = "AVX-512 instructions"
    elif "avx2" in capabilities:
        cpu_flags["avx2"] = "AVX2 instructions"
    elif "avx" in capabilities:
        cpu_flags["avx"] = "AVX instructions"
    elif "sse4_2" in capabilities:
        cpu_flags["sse4_2"] = "SSE4.2 instructions"
    
    if cpu_flags:
        cpu_flags["native"] = "Native CPU (-march=native)"
        cpu_flags["none"] = "No CPU-specific optimizations"
        flag_options["CPU Instructions"] = cpu_flags
    else:
        flag_options["CPU Instructions"] = {
            "native": "Native CPU (-march=native)",
            "none": "No CPU-specific optimizations"
        }
    
    # Advanced optimizations
    advanced_flags = {
        "lto": "Link Time Optimization",
        "pgo": "Profile-Guided Optimization (requires multiple build steps)",
        "none": "No advanced optimizations"
    }
    
    # Add polly only for Clang
    if compiler_type == "clang":
        advanced_flags["polly"] = "LLVM Polly optimizations (Clang only)"
        
    flag_options["Advanced Optimizations"] = advanced_flags
    
    # Debug info
    flag_options["Debug Information"] = {
        "g0": "No debug info",
        "g1": "Minimal debug info",
        "g2": "Standard debug info",
        "g3": "Maximum debug info"
    }
    
    return flag_options

def check_build_tools():
    build_tools = {}
    
    # Check for CMake
    cmake_path = shutil.which("cmake")
    if cmake_path:
        try:
            result = subprocess.run([cmake_path, "--version"],
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.STDOUT,
                                 text=True,
                                 timeout=5)
            version = result.stdout.strip().split('\n')[0]
            build_tools["cmake"] = {"path": cmake_path, "version": version}
            print(f"{Fore.GREEN}Found CMake: {version}{Style.RESET_ALL}")
        except Exception:
            pass
    
    # Check for Ninja
    ninja_path = shutil.which("ninja")
    if ninja_path:
        try:
            result = subprocess.run([ninja_path, "--version"],
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.STDOUT,
                                 text=True,
                                 timeout=5)
            version = result.stdout.strip()
            build_tools["ninja"] = {"path": ninja_path, "version": version}
            print(f"{Fore.GREEN}Found Ninja: {version}{Style.RESET_ALL}")
        except Exception:
            pass
    
    # Check for Make
    make_path = shutil.which("make")
    if make_path:
        try:
            result = subprocess.run([make_path, "--version"],
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.STDOUT,
                                 text=True,
                                 timeout=5)
            version = result.stdout.strip().split('\n')[0]
            build_tools["make"] = {"path": make_path, "version": version}
            print(f"{Fore.GREEN}Found Make: {version}{Style.RESET_ALL}")
        except Exception:
            pass
    
    return build_tools

def main():
    parser = argparse.ArgumentParser(description="DescriptorFactory Build Configuration Tool")
    parser.add_argument("--output", "-o", default="compile_flags.txt", 
                     help="Output file for selected flags (default: compile_flags.txt)")
    parser.add_argument("--verbose", "-v", action="store_true", 
                     help="Show verbose output")
    parser.add_argument("--build-dir", "-b", default="build",
                     help="Build directory name (default: build)")
    parser.add_argument("--skip-deps", "-s", action="store_true",
                     help="Skip dependency checking")
    args = parser.parse_args()
    
    init()
    print(f"{Fore.CYAN}DescriptorFactory Build Configuration Tool{Style.RESET_ALL}")
    print(f"{Fore.YELLOW}Type 'q' at any prompt to exit, '?' for help when available{Style.RESET_ALL}")
    
    # Check for dependencies
    if not args.skip_deps:
        missing_deps = check_dependencies()
        if missing_deps:
            display_install_instructions(missing_deps)
    
    # Detect compiler
    print(f"{Fore.CYAN}Detecting compilers...{Style.RESET_ALL}")
    compilers = get_compiler_info()
    
    if not compilers:
        print(f"{Fore.RED}No suitable compilers found. Please install GCC or Clang.{Style.RESET_ALL}")
        return 1
    
    # Check build tools
    print(f"{Fore.CYAN}Detecting build tools...{Style.RESET_ALL}")
    build_tools = check_build_tools()
    
    if not build_tools:
        print(f"{Fore.RED}No suitable build tools found. Please install CMake, Ninja, or Make.{Style.RESET_ALL}")
        return 1
    
    # Collect compiler choice
    compiler_choices = list(compilers.keys())
    compiler_choice = None
    
    if len(compiler_choices) > 1:
        compiler_choice = get_user_choice("Select compiler", compiler_choices, compiler_choices[0])
    else:
        compiler_choice = compiler_choices[0]
        print(f"{Fore.CYAN}Using compiler: {compiler_choice}{Style.RESET_ALL}")
    
    compiler_type = compilers[compiler_choice]["type"]
    compiler_version = compilers[compiler_choice]["version"]
    print(f"{Fore.GREEN}Selected {compiler_choice}: {compiler_version}{Style.RESET_ALL}")
    
    # Detect CPU capabilities
    print(f"{Fore.CYAN}Detecting CPU capabilities...{Style.RESET_ALL}")
    capabilities = detect_cpu_capabilities()
    if capabilities:
        print(f"{Fore.GREEN}Detected CPU features: {', '.join(capabilities)}{Style.RESET_ALL}")
    else:
        print(f"{Fore.YELLOW}No specific CPU features detected{Style.RESET_ALL}")
    
    # Get available flags
    flag_options = get_optimization_flags(compiler_type, capabilities)
    
    # Collect user choices
    selected_flags = []
    
    # Optimization level
    opt_level = get_user_choice(
        "Select optimization level", 
        list(flag_options["Optimization Level"].keys()),
        "O3",
        help_key="optimization"
    )
    selected_flags.append(f"-{opt_level}")
    
    # CPU instructions
    cpu_choice = get_user_choice(
        "Select CPU instruction set", 
        list(flag_options["CPU Instructions"].keys()),
        "native" if flag_options["CPU Instructions"].get("native") else "none",
        help_key="cpu"
    )
    
    if cpu_choice != "none":
        if cpu_choice == "native":
            selected_flags.append("-march=native")
        elif cpu_choice == "avx512":
            selected_flags.extend(["-mavx512f", "-mavx512dq", "-mavx512bw", "-mavx512vl"])
        elif cpu_choice in ["avx", "avx2", "sse4_2"]:
            selected_flags.append(f"-m{cpu_choice}")
            if cpu_choice == "avx2" and "fma" in capabilities:
                selected_flags.append("-mfma")
    
    # Advanced optimizations
    adv_choice = get_user_choice(
        "Select advanced optimizations", 
        list(flag_options["Advanced Optimizations"].keys()),
        "none",
        help_key="advanced"
    )
    
    if adv_choice != "none":
        if adv_choice == "lto":
            # Use appropriate LTO flag based on compiler
            if compiler_type == "clang":
                selected_flags.append("-flto=thin")
            else:
                selected_flags.append("-flto")
        elif adv_choice == "pgo":
            print(f"{Fore.CYAN}PGO selected - you'll need to run multiple build steps{Style.RESET_ALL}")
            if compiler_type == "clang":
                print("Generate step: -fprofile-instr-generate")
                print("Use step: -fprofile-instr-use=default.profdata")
            else:
                print("Generate step: -fprofile-generate=.")
                print("Use step: -fprofile-use=.")
        elif adv_choice == "polly" and compiler_type == "clang":
            selected_flags.extend(["-mllvm", "-polly"])
    
    # Debug info
    debug_choice = get_user_choice(
        "Select debug information level", 
        list(flag_options["Debug Information"].keys()),
        "g0",
        help_key="debug"
    )
    
    if debug_choice != "g0":
        selected_flags.append(f"-{debug_choice}")
    
    # Additional flags
    print(f"{Fore.CYAN}Common additional flags:{Style.RESET_ALL}")
    
    if get_user_choice("Add -DNDEBUG (disable assert)", ["Yes", "No"], "Yes") == "Yes":
        selected_flags.append("-DNDEBUG")
    
    if get_user_choice("Add -Wall (enable all warnings)", ["Yes", "No"], "Yes") == "Yes":
        selected_flags.append("-Wall")
        selected_flags.append("-Wno-reorder")
    
    if get_user_choice("Add -Wextra (enable extra warnings)", ["Yes", "No"], "No") == "Yes":
        selected_flags.append("-Wextra")
    
    # DescriptorFactory specific options
    print(f"{Fore.CYAN}DescriptorFactory specific options:{Style.RESET_ALL}")
    
    # GUI option
    build_gui = get_user_choice(
        "Build the GUI application", 
        ["Yes", "No"], 
        "Yes",
        help_key="gui"
    ) == "Yes"
    
    # Enable PGO
    enable_pgo = False
    if adv_choice == "pgo":
        enable_pgo = True
        print(f"{Fore.CYAN}Selected PGO - will enable ENABLE_PGO option{Style.RESET_ALL}")
    
    # Build system selection
    build_system = "cmake"
    generator = "Unix Makefiles"
    
    if len(build_tools) > 1:
        build_choices = []
        if "cmake" in build_tools:
            if "ninja" in build_tools:
                build_choices.append("CMake+Ninja")
            if "make" in build_tools:
                build_choices.append("CMake+Make")
        if "make" in build_tools:
            build_choices.append("Make")
        
        build_choice = get_user_choice(
            "Select build system", 
            build_choices, 
            build_choices[0],
            help_key="build_system"
        )
        
        if build_choice == "CMake+Ninja":
            build_system = "cmake"
            generator = "Ninja"
        elif build_choice == "CMake+Make":
            build_system = "cmake"
            generator = "Unix Makefiles"
        elif build_choice == "Make":
            build_system = "make"
    
    # Show final flags
    print(f"{Fore.GREEN}Selected flags:{Style.RESET_ALL}")
    for flag in selected_flags:
        print(f"  {flag}")
    
    # Save to file
    try:
        # Join flags with spaces to avoid newlines that cause CMake issues
        flags_str = " ".join(selected_flags)
        
        with open(args.output, "w") as f:
            f.write(flags_str)
        print(f"{Fore.GREEN}Flags saved to {args.output}{Style.RESET_ALL}")
        
        # Create build script
        build_dir = args.build_dir
        script_path = "build.sh"
        
        cmake_options = []
        if enable_pgo:
            cmake_options.append("-DENABLE_PGO=ON")
        if not build_gui:
            cmake_options.append("-DBUILD_GUI=OFF")
        
        with open(script_path, "w") as f:
            f.write("#!/bin/bash\n\n")
            f.write("# Script generated by DescriptorFactory Build Configuration Tool\n\n")
            
            f.write("# Create build directory\n")
            f.write(f"mkdir -p {build_dir} && cd {build_dir}\n\n")
            
            if build_system == "cmake":
                f.write("# Configure with CMake\n")
                cmake_cmd = f'cmake -G "{generator}" -DCMAKE_CXX_FLAGS="$(cat ../{args.output})" {" ".join(cmake_options)} ..'
                f.write(f"{cmake_cmd}\n\n")
                
                f.write("# Build\n")
                if generator == "Ninja":
                    f.write("ninja\n")
                else:
                    # Use processor count for make
                    f.write('make -j"$(nproc)"\n')
            else:
                f.write("# Build with Make\n")
                f.write(f'CXXFLAGS="$(cat ../{args.output})" make -j"$(nproc)"\n')
            
            # Make script executable
            os.chmod(script_path, 0o755)
        
        print(f"{Fore.GREEN}Build script created: {script_path}{Style.RESET_ALL}")
        
        # Show guidance
        print("\nBuild instructions:")
        print(f"  1. Run the generated script: {Style.BRIGHT}./build.sh{Style.RESET_ALL}")
        print(f"  2. The build will be created in the '{build_dir}' directory")
        
        if enable_pgo:
            print(f"\n{Fore.YELLOW}Profile-Guided Optimization (PGO) workflow:{Style.RESET_ALL}")
            print("  1. Build instrumented version with current settings")
            print("  2. Run benchmarks/tests to generate profile data")
            print("  3. Rebuild using collected profile data")
            
        print(f"\nTo use flags manually:")
        print(f"  • With CMake: cmake -DCMAKE_CXX_FLAGS=\"$(cat {args.output})\" ..")
        print(f"  • With Make: CXXFLAGS=\"$(cat {args.output})\" make")
        
        # Add dependency install commands for reference
        print(f"\n{Fore.CYAN}Dependency installation commands:{Style.RESET_ALL}")
        if sys.platform.startswith("linux"):
            print("Ubuntu/Debian:")
            print("sudo apt install build-essential qtbase5-dev qt5-qmake cmake qtbase5-dev-tools ninja-build cmake libtbb-dev libeigen3-dev librdkit-dev libtbb-dev rapidjson-dev libcxxopts-dev")
        elif sys.platform == "darwin":
            print("macOS:")
            print("brew install cmake eigen rdkit cxxopts tbb qt@5 rapidjson")
        
    except Exception as e:
        print(f"{Fore.RED}Failed to save configuration: {e}{Style.RESET_ALL}")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())