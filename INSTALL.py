#!/usr/bin/env python3
import os
import sys
import argparse
import subprocess
import shutil
import time
import signal

# Simple logger without external dependencies
class Logger:
    # ANSI color codes
    RED = '\033[0;31m'
    GREEN = '\033[0;32m'
    YELLOW = '\033[1;33m'
    BLUE = '\033[0;34m'
    CYAN = '\033[0;36m'
    WHITE = '\033[0;37m'
    NC = '\033[0m'  # No Color
    BOLD = '\033[1m'
    
    def __init__(self, verbose=False):
        self.verbose = verbose
        self.progress_items = 0
        self.progress_current = 0
        self.progress_desc = ""
        
    def _format_message(self, level, message, color):
        return f"{color}[{level}]{self.NC} {message}"
    
    def info(self, message):
        print(self._format_message("INFO", message, self.BLUE))
        
    def debug(self, message):
        if self.verbose:
            print(self._format_message("DEBUG", message, self.WHITE))
        
    def warn(self, message):
        print(self._format_message("WARN", message, self.YELLOW))
        
    def error(self, message):
        print(self._format_message("ERROR", message, self.RED))
        
    def success(self, message):
        print(self._format_message("SUCCESS", message, self.GREEN))
    
    def start_progress(self, total, desc="Processing"):
        self.progress_items = total
        self.progress_current = 0
        self.progress_desc = desc
        print(f"{self.BLUE}[INFO]{self.NC} {desc}: 0/{total} (0.00%)")
        return self
    
    def update(self, n=1):
        self.progress_current += n
        percentage = (self.progress_current / self.progress_items) * 100
        print(f"{self.BLUE}[INFO]{self.NC} {self.progress_desc}: {self.progress_current}/{self.progress_items} ({percentage:.2f}%)")
    
    def set_description(self, desc):
        self.progress_desc = desc
    
    def finish(self):
        if self.progress_items > 0:
            print(f"{self.BLUE}[INFO]{self.NC} {self.progress_desc}: {self.progress_items}/{self.progress_items} (100.00%)")
            self.progress_items = 0
    
    def __enter__(self):
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        self.finish()

def detect_compiler_flags(logger, quick=False):
    """Detect which compiler flags are supported for PGO"""
    logger.info("Detecting supported compiler flags for PGO...")
    
    # If we're in quick mode, just use the most common flags
    if quick:
        logger.info("Quick mode: using standard GCC PGO flags")
        return "-fprofile-generate", "-fprofile-use"
    
    # Create a temporary directory for testing
    temp_dir = "/tmp/compiler_test"
    os.makedirs(temp_dir, exist_ok=True)
    
    # Write a simple test program
    with open(os.path.join(temp_dir, "test.cpp"), "w") as f:
        f.write("#include <iostream>\nint main() { std::cout << \"Hello, World!\" << std::endl; return 0; }")
    
    # Detect compiler and version first
    compiler_info = {}
    for compiler in ["c++", "g++", "clang++"]:
        try:
            result = subprocess.run([compiler, "--version"], 
                                    stdout=subprocess.PIPE, 
                                    stderr=subprocess.PIPE, 
                                    text=True,
                                    timeout=5)
            if result.returncode == 0:
                compiler_info[compiler] = result.stdout.lower()
                logger.info(f"Found compiler: {compiler}")
                logger.debug(f"Version info: {result.stdout.strip()}")
        except subprocess.TimeoutExpired:
            logger.warn(f"Timeout while checking {compiler} version")
        except:
            pass
    
    # Determine which compiler to use
    compiler_cmd = "c++" # default
    if "g++" in compiler_info:
        compiler_cmd = "g++"
    elif "clang++" in compiler_info:
        compiler_cmd = "clang++"
    
    # Get detailed info about chosen compiler
    chosen_version = compiler_info.get(compiler_cmd, "")
    is_clang = "clang" in chosen_version
    is_gcc = "gcc" in chosen_version or "g++" in chosen_version and not is_clang
    
    logger.info(f"Using {compiler_cmd} for PGO flag detection")
    
    # All possible PGO flag combinations to try, in order of preference
    pgo_flag_options = []
    
    # GCC options
    if is_gcc:
        logger.info("Detected GCC-compatible compiler")
        # Try more advanced optimization combinations first - newer GCC with corrections
        pgo_flag_options.extend([
            # Standard with corrections
            ("-fprofile-generate", "-fprofile-use -fprofile-correction"),
            # Enhanced with directory specification
            ("-fprofile-generate=.", "-fprofile-use=. -fprofile-correction"),
            # With branch probabilities
            ("-fprofile-generate", "-fprofile-use -fbranch-probabilities"),
            # With corrections and probabilities
            ("-fprofile-generate", "-fprofile-use -fprofile-correction -fbranch-probabilities"),
            # Standard GCC PGO
            ("-fprofile-generate", "-fprofile-use"),
            # GCC with path
            ("-fprofile-generate=.", "-fprofile-use=."),
        ])
    
    # Clang options
    if is_clang:
        logger.info("Detected Clang-compatible compiler")
        pgo_flag_options.extend([
            # Standard Clang LLVM profiling
            ("-fprofile-instr-generate", "-fprofile-instr-use=default.profdata"),
            # With merging
            ("-fprofile-instr-generate", "-fprofile-instr-use=default.profdata -Wno-profile-instr-out-of-date"),
            # With specified file
            ("-fprofile-instr-generate=profile.profraw", "-fprofile-instr-use=profile.profdata"),
            # Basic Clang instrumentation
            ("-fprofile-instr-generate", "-fprofile-instr-use"),
            # Old style GCC compatible
            ("-fprofile-generate", "-fprofile-use"),
        ])
    
    # If no compiler was detected, add default options
    if not pgo_flag_options:
        logger.warn("Could not identify compiler type, trying generic options")
        pgo_flag_options.extend([
            ("-fprofile-generate", "-fprofile-use"),
            ("-fprofile-instr-generate", "-fprofile-instr-use"),
        ])
    
    # Test each flag combination
    gen_flag, use_flag = "", ""
    for gen, use in pgo_flag_options:
        # Try to compile with generation flag
        compile_cmd = f"{compiler_cmd} {gen} -o test test.cpp"
        logger.debug(f"Testing compilation with: {compile_cmd}")
        
        try:
            compile_result = subprocess.run(
                compile_cmd, 
                shell=True,
                cwd=temp_dir,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                timeout=10
            )
            
            if compile_result.returncode == 0:
                logger.debug(f"Compilation with {gen} succeeded")
                
                # Run the compiled binary to generate profile data
                run_result = subprocess.run(
                    "./test", 
                    shell=True,
                    cwd=temp_dir,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True,
                    timeout=5
                )
                
                # Now try to use the profile data
                if run_result.returncode == 0:
                    use_cmd = f"{compiler_cmd} {use} -o test_optimized test.cpp"
                    logger.debug(f"Testing optimization with: {use_cmd}")
                    
                    use_result = subprocess.run(
                        use_cmd,
                        shell=True,
                        cwd=temp_dir,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE,
                        text=True,
                        timeout=10
                    )
                    
                    if use_result.returncode == 0:
                        gen_flag, use_flag = gen, use
                        logger.success(f"Found working PGO flags: {gen} → {use}")
                        break
                
        except subprocess.TimeoutExpired:
            logger.warn(f"Timeout while testing {gen}")
        except Exception as e:
            logger.debug(f"Error testing {gen}: {str(e)}")
    
    # Clean up
    try:
        shutil.rmtree(temp_dir, ignore_errors=True)
    except:
        logger.debug(f"Could not clean up {temp_dir}")
    
    # If we couldn't find a working solution, fallback to CMake's own mechanism
    if not gen_flag:
        logger.warn("Could not determine working PGO flags, falling back to CMake defaults")
        return None, None
    
    return gen_flag, use_flag

def run_command(cmd, logger, cwd=None, env=None, show_output=False, shell=False, timeout=None, hide_warnings=True):
    logger.debug(f"Running command: {cmd if shell else ' '.join(cmd)}")
    try:
        process = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            cwd=cwd,
            env=env,
            shell=shell
        )
        
        stdout_data = []
        stderr_data = []
        
        # Non-blocking reads with timeout
        import select
        if timeout:
            end_time = time.time() + timeout
        
        # Set up polling for stdout and stderr
        poller = select.poll()
        poller.register(process.stdout, select.POLLIN)
        poller.register(process.stderr, select.POLLIN)
        
        # Handle output with polling to avoid blocking
        while process.poll() is None:
            try:
                # Check for timeout
                if timeout and time.time() > end_time:
                    process.terminate()
                    raise subprocess.TimeoutExpired(cmd, timeout)
                
                # Poll with a short timeout
                events = poller.poll(100)  # 100ms timeout
                
                for fd, event in events:
                    if fd == process.stdout.fileno():
                        line = process.stdout.readline()
                        if line:
                            # Save all output
                            stdout_data.append(line)
                            # But only display if it's not a warning or if we're showing warnings
                            if show_output or logger.verbose:
                                if not hide_warnings or "warning:" not in line:
                                    print(line.rstrip())
                    elif fd == process.stderr.fileno():
                        line = process.stderr.readline()
                        if line:
                            # Save all output
                            stderr_data.append(line)
                            # But only display if it's not a warning or if we're showing warnings
                            if show_output or logger.verbose:
                                if not hide_warnings or "warning:" not in line:
                                    print(line.rstrip())
            except KeyboardInterrupt:
                # Handle Ctrl+C gracefully
                process.terminate()
                logger.warn("Command interrupted by user")
                return False, "".join(stdout_data), "Command interrupted by user"
            
            # Sleep briefly to avoid high CPU usage
            time.sleep(0.001)
        
        # Collect any remaining output
        for line in process.stdout:
            stdout_data.append(line)
            if show_output or logger.verbose:
                if not hide_warnings or "warning:" not in line:
                    print(line.rstrip())
            
        for line in process.stderr:
            stderr_data.append(line)
            if show_output or logger.verbose:
                if not hide_warnings or "warning:" not in line:
                    print(line.rstrip())
        
        stdout = "".join(stdout_data)
        stderr = "".join(stderr_data)
        
        if process.returncode != 0:
            logger.error(f"Command failed with exit code {process.returncode}")
            if not show_output and not logger.verbose:
                logger.debug(f"Command output:\n{stdout}")
                logger.error(f"Command error:\n{stderr}")
            return False, stdout, stderr
        
        return True, stdout, stderr
    except subprocess.TimeoutExpired:
        logger.error(f"Command timed out after {timeout} seconds")
        return False, "", "Command timed out"
    except Exception as e:
        logger.error(f"Exception when running command: {e}")
        return False, "", str(e)

def ensure_input_directory(input_dir, logger):
    if not os.path.exists(input_dir):
        logger.info(f"Creating input directory: {input_dir}")
        os.makedirs(input_dir, exist_ok=True)
    
    # Check if directory is empty or has no CSV files
    if not any(f.endswith('.csv') for f in os.listdir(input_dir)):
        logger.info("No CSV files found in input directory. Creating sample data file...")
        
        # Use proper SMILES header (uppercase) to match program expectations
        sample_data = """SMILES
CC(C)CC1=CC=C(C=C1)C(C)C(=O)O
CC(=O)OC1=CC=CC=C1C(=O)O
CC1=C(C=C(C=C1)CC(=O)O)OC
CC1=CC(=C(C(=C1)OC)O)O
CC1=CC=C(O1)C(=O)O
CC(C)C1=CC=CC(=C1)C(=O)O
CC(C)(C)C1=CC(=C(C=C1)O)CC(=O)O
CC(C)(C)C1=CC(=CC(=C1)O)CC(=O)O
CC(C)C1=CC=CC=C1CC(=O)O
CC(C)C1=CC=CC=C1CCC(=O)O"""
        
        sample_file_path = os.path.join(input_dir, "sample_smiles.csv")
        with open(sample_file_path, 'w') as f:
            f.write(sample_data)
        
        logger.success(f"Created sample SMILES file at {sample_file_path}")

def stage_instrumented_build(build_dir, jobs, logger, portable=True, pgo_flags=None, timeout=None, hide_warnings=True):
    logger.info("STAGE 1: Building instrumented binary")
    
    # Create or clean build directory
    if not os.path.exists(build_dir):
        os.makedirs(build_dir)
    
    # Clean any existing build artifacts
    for item in os.listdir(build_dir):
        item_path = os.path.join(build_dir, item)
        if os.path.isfile(item_path):
            os.unlink(item_path)
        elif os.path.isdir(item_path):
            shutil.rmtree(item_path)
    
    # Configure with instrumentation
    cmake_command = [
        "cmake", 
        "-GNinja", 
        "-DCMAKE_BUILD_TYPE=Release",
        "-DENABLE_PGO=ON",
        "-DGENERATE_PROFILE=ON",
        "-DUSE_PROFILE=OFF",
        # Disable warnings we don't care about
        "-DCMAKE_CXX_FLAGS=-Wno-unused-variable -Wno-unused-function -Wno-unused-parameter -Wno-sign-compare -Wno-reorder"
    ]
    
    # Add custom PGO flags if detected
    if pgo_flags and pgo_flags[0]:
        # Append PGO flags to existing flags
        cmake_command[-1] = f"{cmake_command[-1]} {pgo_flags[0]}"
    
    # Use portable build by default to avoid SIGILL
    if portable:
        cmake_command.append("-DOPTIMIZE_FOR_NATIVE=OFF")
    
    cmake_command.append("..")
    
    success, _, _ = run_command(cmake_command, logger, cwd=build_dir, show_output=True, timeout=timeout, hide_warnings=hide_warnings)
    if not success:
        logger.error("Failed to configure instrumented build")
        return False
    
    # Build instrumented binary
    build_command = ["ninja", f"-j{jobs}"]
    success, _, _ = run_command(build_command, logger, cwd=build_dir, show_output=True, timeout=timeout, hide_warnings=hide_warnings)
    if not success:
        logger.error("Failed to build instrumented binary")
        return False
    
    # Verify binary exists
    binary_path = os.path.join(build_dir, "desfact")
    if not os.path.exists(binary_path):
        logger.error(f"Binary not found at {binary_path} after build")
        return False
    
    logger.success(f"Successfully built instrumented binary at {binary_path}")
    return True

def stage_generate_profile(build_dir, input_dir, logger, timeout=None, hide_warnings=True):
    logger.info("STAGE 2: Running benchmarks to generate profile data")
    
    # Check if instrumented binary exists
    binary_path = os.path.join(build_dir, "desfact")
    if not os.path.exists(binary_path):
        logger.error(f"Instrumented binary not found at {binary_path}")
        return False
    
    # Check if binary is executable
    if not os.access(binary_path, os.X_OK):
        logger.warn(f"Binary at {binary_path} is not executable, attempting to set permissions")
        try:
            os.chmod(binary_path, 0o755)  # rwxr-xr-x
        except Exception as e:
            logger.error(f"Failed to set executable permissions: {e}")
            return False
    
    # Find all CSV files in the input directory
    csv_files = [f for f in os.listdir(input_dir) if f.endswith('.csv')]
    if not csv_files:
        logger.error(f"No CSV files found in {input_dir}")
        return False
    
    logger.info(f"Found {len(csv_files)} CSV files for benchmark")
    
    # Before running benchmarks, add:
    os.makedirs(os.path.join(build_dir, "output"), exist_ok=True)
    
    # Run benchmark for each CSV file
    with logger.start_progress(len(csv_files), "Running benchmarks") as progress:
        for i, csv_file in enumerate(csv_files):
            progress.set_description(f"Benchmark {i+1}/{len(csv_files)}")
            csv_path = os.path.join(os.path.abspath(input_dir), csv_file)
            output_path = os.path.abspath(os.path.join(build_dir, "output", f"temp_output_{i}.csv"))
            
            # Properly format the command with all required arguments
            benchmark_command = [
                "./desfact",
                "-i", csv_path,
                "-o", output_path,
                "-d", "all"  # Use all descriptors for better profiling
            ]
            
            logger.debug(f"Running benchmark: {' '.join(benchmark_command)}")
            
            # Run the benchmark
            success, stdout, stderr = run_command(benchmark_command, logger, cwd=build_dir, timeout=timeout, hide_warnings=hide_warnings)
            if not success:
                logger.warn(f"Benchmark failed with default arguments: {stderr}")
                
                # If that fails, try specifying the SMILES column name explicitly
                logger.warn("Retrying with explicit SMILES column name")
                benchmark_command_alt = [
                    "./desfact",
                    "-i", csv_path,
                    "-o", output_path,
                    "-d", "all",
                    "-s", "SMILES",  # Explicitly specify SMILES column name
                    "--verbose"      # Add verbose output for debugging
                ]
                success, stdout, stderr = run_command(benchmark_command_alt, logger, cwd=build_dir, timeout=timeout, hide_warnings=hide_warnings)
                if not success:
                    logger.error(f"Benchmark also failed with explicit SMILES column: {stderr}")
                    
                    # Try with simpler descriptor set
                    logger.warn("Trying with minimal descriptor set")
                    simple_cmd = [
                        "./desfact",
                        "-i", csv_path,
                        "-o", output_path,
                        "-d", "all"
                    ]
                    success, stdout, stderr = run_command(simple_cmd, logger, cwd=build_dir, timeout=timeout, hide_warnings=hide_warnings)
                    if not success:
                        logger.error(f"All benchmark attempts failed for {csv_file}")
                        continue
            
            # Clean up temporary output file if it exists
            try:
                if os.path.exists(output_path):
                    os.unlink(output_path)
            except:
                pass
                
            progress.update()
    
    # Look for profile data (recursive search)
    profile_files = []
    for root, dirs, files in os.walk(build_dir):
        for file in files:
            if file.endswith('.gcda') or file.endswith('.profraw'):
                profile_files.append(os.path.join(root, file))
    
    if not profile_files:
        logger.error("No profile data (.gcda or .profraw) was generated")
        return False
    
    # Make sure we have at least some profile data
    logger.success(f"Successfully generated profile data ({len(profile_files)} files)")
    
    # Write profile file list for debugging
    with open(os.path.join(build_dir, "profile_files.txt"), "w") as f:
        for profile_file in profile_files:
            f.write(f"{profile_file}\n")
    
    return True

def stage_optimized_build(build_dir, jobs, logger, portable=True, pgo_flags=None, timeout=None, hide_warnings=True):
    logger.info("STAGE 3: Building optimized binary with profile data")
    
    # Ensure profile data is available
    profile_files = []
    for root, dirs, files in os.walk(build_dir):
        for file in files:
            if file.endswith('.gcda') or file.endswith('.profraw'):
                profile_files.append(os.path.join(root, file))
    
    if not profile_files:
        logger.error("No profile data found. Cannot build optimized binary.")
        return False
    
    logger.info(f"Found {len(profile_files)} profile data files")
    
    # Configure with profile use
    cmake_command = [
        "cmake", 
        "-GNinja",
        "-DCMAKE_BUILD_TYPE=Release",
        "-DENABLE_PGO=ON",
        "-DGENERATE_PROFILE=OFF",
        "-DUSE_PROFILE=ON",
        # Disable warnings we don't care about
        "-DCMAKE_CXX_FLAGS=-Wno-unused-variable -Wno-unused-function -Wno-unused-parameter -Wno-sign-compare -Wno-reorder"
    ]
    
    # Special handling for GCC PGO - create an empty merged profile if needed
    try:
        # Check for GCC style profile data
        if any(f.endswith('.gcda') for f in profile_files):
            logger.info("Found GCC-style profile data (.gcda files)")
            
            # Add additional options for PGO
            cmake_command.append("-DCMAKE_EXE_LINKER_FLAGS=-fprofile-use -fprofile-correction")
            
            # Make sure profile data is recognized
            if pgo_flags and pgo_flags[1]:
                # Append PGO flags to existing flags with explicit directory
                cmake_command[-2] = f"{cmake_command[-2]} {pgo_flags[1]}"
        
        # If Clang-style, we need to merge profiles
        elif any(f.endswith('.profraw') for f in profile_files):
            logger.info("Found Clang-style profile data (.profraw files)")
            # For Clang we might need to merge the profile data
            merge_cmd = ["llvm-profdata", "merge", "-output=default.profdata"]
            merge_cmd.extend([f for f in profile_files if f.endswith('.profraw')])
            
            run_command(merge_cmd, logger, cwd=build_dir)
            logger.info("Merged profile data files")
    except Exception as e:
        logger.warn(f"Error preparing profile data: {e}")
    
    # Add custom PGO flags if detected
    if pgo_flags and pgo_flags[1]:
        # Append PGO flags to existing flags
        cmake_command[-1] = f"{cmake_command[-1]} {pgo_flags[1]}"
    
    # Use portable build by default to avoid SIGILL
    if portable:
        cmake_command.append("-DOPTIMIZE_FOR_NATIVE=OFF")
    
    cmake_command.append("..")
    
    success, _, _ = run_command(cmake_command, logger, cwd=build_dir, show_output=True, timeout=timeout, hide_warnings=hide_warnings)
    if not success:
        logger.error("Failed to configure optimized build")
        return False
    
    # Build optimized binary
    build_command = ["ninja", f"-j{jobs}"]
    success, _, _ = run_command(build_command, logger, cwd=build_dir, show_output=True, timeout=timeout, hide_warnings=hide_warnings)
    if not success:
        logger.error("Failed to build optimized binary")
        return False
    
    # Verify binary exists
    binary_path = os.path.join(build_dir, "desfact")
    if not os.path.exists(binary_path):
        logger.error(f"Binary not found at {binary_path} after build")
        return False
    
    logger.success("Successfully built optimized binary")
    return True

def signal_handler(sig, frame):
    print("\nInterrupted by user. Exiting gracefully...")
    sys.exit(1)

def main():
    # Register signal handler for clean exits
    signal.signal(signal.SIGINT, signal_handler)
    
    # Fast help display that doesn't import heavy dependencies
    if len(sys.argv) == 2 and sys.argv[1] in ['--help', '-h']:
        print("Profile-Guided Optimization (PGO) Tool for DescriptorFactory")
        print("")
        print("This script automates the Profile-Guided Optimization process:")
        print("  1. Builds an instrumented binary for profile generation")
        print("  2. Runs a benchmark to generate profile data")
        print("  3. Builds an optimized binary using the profile data")
        print("")
        print("Usage:")
        print("  ./run_pgo.py [options]")
        print("")
        print("Options:")
        print("  -h, --help                Display this help message")
        print("  -j, --jobs NUM            Set number of parallel jobs (default: number of CPU cores)")
        print("  -i, --input-dir DIR       Directory with test SMILES files (default: ./input)")
        print("  -b, --build-dir DIR       Build directory (default: ./build)")
        print("  -v, --verbose             Enable verbose output")
        print("  -q, --quick               Quick mode: skip thorough testing of compiler flags")
        print("  -t, --timeout SEC         Timeout for build processes in seconds (default: none)")
        print("  -w, --show-warnings       Show compiler warnings (default: hidden)")
        print("  --native                  Optimize for native CPU (may cause SIGILL on different CPUs)")
        print("  --skip-instrumented       Skip instrumented build stage")
        print("  --skip-profile            Skip profile generation stage")
        print("  --skip-optimized          Skip optimized build stage")
        return 0
    
    parser = argparse.ArgumentParser(
        description="Run Profile-Guided Optimization (PGO) for DescriptorFactory",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument("-j", "--jobs", type=int, default=os.cpu_count(),
                        help="Number of parallel build jobs")
    parser.add_argument("-i", "--input-dir", default="./input",
                        help="Directory containing input CSV files for benchmarking")
    parser.add_argument("-b", "--build-dir", default="./build",
                        help="Build directory")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Enable verbose output")
    parser.add_argument("-q", "--quick", action="store_true",
                        help="Quick mode: skip thorough testing of compiler flags")
    parser.add_argument("-t", "--timeout", type=int, default=None,
                        help="Timeout for build processes in seconds")
    parser.add_argument("-w", "--show-warnings", action="store_true",
                        help="Show compiler warnings (default: hidden)")
    parser.add_argument("--native", action="store_true",
                        help="Optimize for native CPU (may cause SIGILL on different CPUs)")
    parser.add_argument("--skip-instrumented", action="store_true",
                        help="Skip instrumented build stage")
    parser.add_argument("--skip-profile", action="store_true",
                        help="Skip profile generation stage")
    parser.add_argument("--skip-optimized", action="store_true",
                        help="Skip optimized build stage")
    
    try:
        args = parser.parse_args()
    except SystemExit:
        return 1
    
    logger = Logger(verbose=args.verbose)
    
    # Use portable build by default (disable -march=native)
    portable = not args.native
    # Whether to hide compiler warnings
    hide_warnings = not args.show_warnings
    
    logger.info("Starting Profile-Guided Optimization process")
    logger.debug(f"Jobs: {args.jobs}")
    logger.debug(f"Input directory: {args.input_dir}")
    logger.debug(f"Build directory: {args.build_dir}")
    logger.debug(f"Portable build: {portable}")
    logger.debug(f"Quick mode: {args.quick}")
    logger.debug(f"Hide warnings: {hide_warnings}")
    if args.timeout:
        logger.debug(f"Command timeout: {args.timeout}s")
    
    try:
        # Detect compiler flags
        pgo_flags = detect_compiler_flags(logger, args.quick)
        
        # Ensure input directory exists and has data
        ensure_input_directory(args.input_dir, logger)
        
        start_time = time.time()
        success = True
        
        # Stage 1: Build instrumented binary
        if not args.skip_instrumented:
            success = stage_instrumented_build(args.build_dir, args.jobs, logger, portable, 
                                              pgo_flags, args.timeout, hide_warnings)
            if not success:
                logger.error("Failed to build instrumented binary. Aborting.")
                return 1
        else:
            logger.info("Skipping instrumented build stage")
        
        # Stage 2: Generate profile
        if not args.skip_profile and success:
            success = stage_generate_profile(args.build_dir, args.input_dir, logger, 
                                            args.timeout, hide_warnings)
            if not success:
                logger.error("Failed to generate profile data. Aborting.")
                return 1
        else:
            logger.info("Skipping profile generation stage")
        
        # Stage 3: Build optimized binary
        if not args.skip_optimized and success:
            success = stage_optimized_build(args.build_dir, args.jobs, logger, portable, 
                                           pgo_flags, args.timeout, hide_warnings)
            if not success:
                logger.error("Failed to build optimized binary.")
                return 1
        else:
            logger.info("Skipping optimized build stage")
        
        elapsed_time = time.time() - start_time
        
        if success:
            logger.success(f"PGO process completed successfully in {elapsed_time:.2f} seconds")
            logger.info(f"Optimized binary available at: {os.path.join(args.build_dir, 'desfact')}")
            return 0
        else:
            logger.error(f"PGO process failed after {elapsed_time:.2f} seconds")
            return 1
    
    except KeyboardInterrupt:
        logger.warn("Process interrupted by user. Exiting gracefully...")
        return 1
    except Exception as e:
        logger.error(f"Unexpected error: {str(e)}")
        return 1

if __name__ == "__main__":
    sys.exit(main()) 