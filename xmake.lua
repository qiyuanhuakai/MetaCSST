set_project("MetaCSST")
set_version("1.0.1")
set_xmakever("2.8.7")

add_rules("mode.debug", "mode.release")
set_languages("c++20")
set_defaultmode("release")
set_warnings("all")



add_requires("nlohmann_json 3.11.3")
add_requires("toml++ 3.4.0")
add_requires("yaml-cpp 0.8.0")

local src_dir = "src"
local example_dir = "test_baseline/example"
local test_dir = "test_output"

local function compare_cmd(test_file, expected_file, label)
    return string.format(
        "import sys; from collections import Counter; "
            .. "test_lines=Counter(open('%s','r').readlines()); "
            .. "expected_lines=Counter(open('%s','r').readlines()); "
            .. "common=sum((test_lines & expected_lines).values()); "
            .. "total=sum(expected_lines.values()); "
            .. "print(f'%s: {common}/{total} lines ({100*common/total:.1f}%%) [strict line-content multiset equality, order-insensitive]'); "
            .. "sys.exit(0 if test_lines==expected_lines else 1)",
        test_file,
        expected_file,
        label
    )
end

local function consistency_cmd(file_a, file_b, label)
    return string.format(
        "import sys; from collections import Counter; "
            .. "a=Counter(open('%s','r').readlines()); "
            .. "b=Counter(open('%s','r').readlines()); "
            .. "common=sum((a & b).values()); "
            .. "total=max(sum(a.values()), sum(b.values())); "
            .. "print(f'%s: {common}/{total} lines ({100*common/total if total else 100:.1f}%%) [strict line-content multiset equality, order-insensitive]'); "
            .. "sys.exit(0 if a==b else 1)",
        file_a,
        file_b,
        label
    )
end

target("MetaCSSTmain")
    set_kind("binary")
    set_targetdir(".")
    set_basename("MetaCSSTmain")
    add_files("src/main_modern.cpp")
    add_includedirs(src_dir)
    add_packages("nlohmann_json", "toml++", "yaml-cpp")
    add_cxflags("-O2", "-Wall", "-Wextra", "-pthread")
    if is_plat("linux") then
        add_syslinks("pthread", "z")
    end

target("MetaCSSTsub")
    set_kind("binary")
    set_targetdir(".")
    set_basename("MetaCSSTsub")
    add_files("src/sub_modern.cpp")
    add_includedirs(src_dir)
    add_packages("nlohmann_json", "toml++", "yaml-cpp")
    add_cxflags("-O2", "-Wall", "-Wextra", "-pthread")
    if is_plat("linux") then
        add_syslinks("pthread", "z")
    end

task("deps")
    set_menu({
        usage = "xmake deps",
        description = "Resolve parser dependencies using xmake package manager"
    })
    on_run(function ()
        os.execv("xmake", {"f", "-y", "-m", "release"})
    end)

task("modern")
    set_menu({
        usage = "xmake modern",
        description = "Build MetaCSSTmain and MetaCSSTsub natively with xmake"
    })
    on_run(function ()
        os.execv("xmake", {"f", "-y", "-m", "release"})
        os.execv("xmake", {"build", "-a"})
    end)

task("verify-json")
    set_menu({
        usage = "xmake verify-json",
        description = "Run JSON pipeline verification"
    })
    on_run(function ()
        os.execv("xmake", {"modern"})
        os.execv("mkdir", {"-p", test_dir})
        os.execv("./MetaCSSTmain", {"-build", "config.json", "-in", path.join(example_dir, "hv29.fa"), "-out", path.join(test_dir, "raw"), "-thread", "4"})
        os.execv("python3", {path.join(src_dir, "call_vr.py"), path.join(test_dir, "raw", "raw.gtf"), path.join(example_dir, "hv29.fa"), path.join(test_dir, "final.gtf")})
        os.execv("python3", {"-c", compare_cmd(path.join(test_dir, "final.gtf"), path.join(example_dir, "out-DGR.gtf"), "JSON Match")})
    end)

task("verify-toml")
    set_menu({
        usage = "xmake verify-toml",
        description = "Run TOML pipeline verification"
    })
    on_run(function ()
        os.execv("xmake", {"modern"})
        os.execv("mkdir", {"-p", test_dir})
        os.execv("./MetaCSSTmain", {"-build", "config.toml", "-in", path.join(example_dir, "hv29.fa"), "-out", path.join(test_dir, "toml_raw"), "-thread", "1"})
        os.execv("python3", {path.join(src_dir, "call_vr.py"), path.join(test_dir, "toml_raw", "raw.gtf"), path.join(example_dir, "hv29.fa"), path.join(test_dir, "toml_final.gtf")})
        os.execv("python3", {"-c", compare_cmd(path.join(test_dir, "toml_final.gtf"), path.join(example_dir, "out-DGR.gtf"), "TOML Match")})
    end)

task("verify-yaml")
    set_menu({
        usage = "xmake verify-yaml",
        description = "Run YAML pipeline verification"
    })
    on_run(function ()
        os.execv("xmake", {"modern"})
        os.execv("mkdir", {"-p", test_dir})
        os.execv("./MetaCSSTmain", {"-build", "config.yaml", "-in", path.join(example_dir, "hv29.fa"), "-out", path.join(test_dir, "yaml_raw"), "-thread", "1"})
        os.execv("python3", {path.join(src_dir, "call_vr.py"), path.join(test_dir, "yaml_raw", "raw.gtf"), path.join(example_dir, "hv29.fa"), path.join(test_dir, "yaml_final.gtf")})
        os.execv("python3", {"-c", compare_cmd(path.join(test_dir, "yaml_final.gtf"), path.join(example_dir, "out-DGR.gtf"), "YAML Match")})
    end)

task("verify-compressed")
    set_menu({
        usage = "xmake verify-compressed",
        description = "Run compressed FASTA pipeline verification"
    })
    on_run(function ()
        os.execv("xmake", {"modern"})
        os.execv("mkdir", {"-p", test_dir})
        os.execv("./MetaCSSTmain", {"-build", "config.json", "-in", path.join(example_dir, "hv29.fa.gz"), "-out", path.join(test_dir, "gz_raw"), "-thread", "1"})
        os.execv("python3", {path.join(src_dir, "call_vr.py"), path.join(test_dir, "gz_raw", "raw.gtf"), path.join(example_dir, "hv29.fa"), path.join(test_dir, "gz_final.gtf")})
        os.execv("python3", {"-c", compare_cmd(path.join(test_dir, "gz_final.gtf"), path.join(example_dir, "out-DGR.gtf"), "GZ Match")})
    end)

task("verify-thread-consistency")
    set_menu({
        usage = "xmake verify-thread-consistency",
        description = "Verify main pipeline thread consistency"
    })
    on_run(function ()
        os.execv("xmake", {"modern"})
        os.execv("mkdir", {"-p", test_dir})
        os.execv("./MetaCSSTmain", {"-build", "config.json", "-in", path.join(example_dir, "hv29.fa"), "-out", path.join(test_dir, "raw_t1"), "-thread", "1"})
        os.execv("./MetaCSSTmain", {"-build", "config.json", "-in", path.join(example_dir, "hv29.fa"), "-out", path.join(test_dir, "raw_t4"), "-thread", "4"})
        os.execv("python3", {"-c", consistency_cmd(path.join(test_dir, "raw_t1", "raw.gtf"), path.join(test_dir, "raw_t4", "raw.gtf"), "MAIN Thread Consistency")})
    end)

task("verify-sub-consistency")
    set_menu({
        usage = "xmake verify-sub-consistency",
        description = "Verify sub pipeline thread consistency"
    })
    on_run(function ()
        os.execv("xmake", {"modern"})
        os.execv("mkdir", {"-p", test_dir})
        os.execv("./MetaCSSTsub", {"-build", "config.json", "-in", path.join(example_dir, "hv29.fa"), "-out", path.join(test_dir, "sub_t1"), "-thread", "1"})
        os.execv("./MetaCSSTsub", {"-build", "config.json", "-in", path.join(example_dir, "hv29.fa"), "-out", path.join(test_dir, "sub_t2"), "-thread", "2"})
        os.execv("python3", {"-c", consistency_cmd(path.join(test_dir, "sub_t1", "out.txt"), path.join(test_dir, "sub_t2", "out.txt"), "SUB Thread Consistency")})
    end)

task("verify")
    set_menu({
        usage = "xmake verify",
        description = "Run all integration verification pipelines (JSON, TOML, YAML, compressed, thread consistency)"
    })
    on_run(function ()
        os.execv("xmake", {"verify-json"})
        os.execv("xmake", {"verify-toml"})
        os.execv("xmake", {"verify-yaml"})
        os.execv("xmake", {"verify-compressed"})
        os.execv("xmake", {"verify-thread-consistency"})
        os.execv("xmake", {"verify-sub-consistency"})
        print("All integration verify pipelines passed.")
    end)

target("test_config")
    set_kind("binary")
    set_targetdir("test")
    set_basename("test_config")
    add_files("test/test_config.cpp")
    add_includedirs(src_dir)
    add_packages("nlohmann_json", "toml++", "yaml-cpp")
    add_cxflags("-O2", "-Wall", "-Wextra", "-pthread")
    if is_plat("linux") then
        add_syslinks("pthread", "z")
    end

target("test_fun")
    set_kind("binary")
    set_targetdir("test")
    set_basename("test_fun")
    add_files("test/test_fun.cpp")
    add_includedirs(src_dir)
    add_packages("nlohmann_json", "toml++", "yaml-cpp")
    add_cxflags("-O2", "-Wall", "-Wextra", "-pthread")
    if is_plat("linux") then
        add_syslinks("pthread", "z")
    end

target("test_fasta")
    set_kind("binary")
    set_targetdir("test")
    set_basename("test_fasta")
    add_files("test/test_fasta.cpp")
    add_includedirs(src_dir)
    add_packages("nlohmann_json", "toml++", "yaml-cpp")
    add_cxflags("-O2", "-Wall", "-Wextra", "-pthread")
    if is_plat("linux") then
        add_syslinks("pthread", "z")
    end

target("test_common")
    set_kind("binary")
    set_targetdir("test")
    set_basename("test_common")
    add_files("test/test_common.cpp")
    add_includedirs(src_dir)
    add_packages("nlohmann_json", "toml++", "yaml-cpp")
    add_cxflags("-O2", "-Wall", "-Wextra", "-pthread")
    if is_plat("linux") then
        add_syslinks("pthread", "z")
    end

task("test")
    set_menu({
        usage = "xmake test",
        description = "Build and run all unit tests from test/ directory"
    })
    on_run(function ()
        os.execv("xmake", {"build", "test_config"})
        os.execv("xmake", {"build", "test_fun"})
        os.execv("xmake", {"build", "test_fasta"})
        os.execv("xmake", {"build", "test_common"})

        print("\n========================================")
        print("Running Unit Tests")
        print("========================================\n")

        local test_dir = "test"
        local has_failure = false
        
        local tests = {"test_config", "test_fun", "test_fasta", "test_common"}
        for _, test_name in ipairs(tests) do
            print("\n>>> Running " .. test_name .. "...")
            local ok = try {
                function ()
                    os.execv(path.join(test_dir, test_name))
                    return true
                end,
                catch {
                    function (errors)
                        return false
                    end
                }
            }
            if not ok then
                has_failure = true
                print("FAILED: " .. test_name)
            end
        end
        
        if has_failure then
            print("\n========================================")
            print("SOME UNIT TESTS FAILED")
            print("========================================")
            os.exit(1)
        else
            print("\n========================================")
            print("All unit tests passed!")
            print("========================================")
        end
    end)


task("example")
    set_menu({
        usage = "xmake example",
        description = "Print runnable MetaCSST usage examples"
    })
    on_run(function ()
        print("MetaCSST usage examples:")
        print("  1) Build models + run full prediction (JSON config)")
        print("     ./MetaCSSTmain -build config.json -in example/hv29.fa -out example/run_out -thread 4")
        print("  2) Build models + run with TOML config")
        print("     ./MetaCSSTmain -build config.toml -in example/hv29.fa -out example/run_out_toml -thread 4")
        print("  3) Build models + run with YAML config")
        print("     ./MetaCSSTmain -build config.yaml -in example/hv29.fa -out example/run_out_yaml -thread 4")
        print("  4) Sub-structure scan only (TR/VR/RT)")
        print("     ./MetaCSSTsub -build config.json -in example/hv29.fa -out example/sub_out -thread 2")
        print("  5) Post-process raw output")
        print("     python3 src/call_vr.py example/run_out/raw.gtf example/hv29.fa example/final.gtf")
    end)

task("pack")
    set_menu({
        usage = "xmake pack",
        description = "Build and package MetaCSST for distribution (handles config->align dependency chain)"
    })
    on_run(function ()
        local dist_version = "1.0.1"
        local dist_name = string.format("MetaCSST-%s-%s-%s", dist_version, os.host(), os.arch())
        local dist_dir = path.join("dist", dist_name)
        local bin_dir = path.join(dist_dir, "bin")
        local config_dir = path.join(dist_dir, "config")
        local align_dir = path.join(config_dir, "align")
        local example_out_dir = path.join(dist_dir, "example")

        os.execv("xmake", {"modern"})

        os.execv("rm", {"-rf", dist_dir})
        os.execv("mkdir", {"-p", bin_dir, config_dir, align_dir, example_out_dir})

        os.execv("cp", {"MetaCSSTmain", bin_dir})
        os.execv("cp", {"MetaCSSTsub", bin_dir})
        os.execv("cp", {path.join(src_dir, "call_vr.py"), bin_dir})

        os.execv("cp", {"config.json", config_dir})
        os.execv("cp", {"config.toml", config_dir})
        os.execv("cp", {"config.yaml", config_dir})

        for _, file in ipairs(os.files("align/*.align")) do
            os.execv("cp", {file, align_dir})
        end

        os.execv("cp", {"README", dist_dir})
        os.execv("cp", {path.join(example_dir, "hv29.fa"), example_out_dir})
        os.execv("cp", {path.join(example_dir, "out-DGR.gtf"), example_out_dir})

        io.writefile(path.join(dist_dir, "requirements.txt"), "numpy\nnumba\npyfastx\n")

        local tarball = dist_dir .. ".tar.gz"
        os.execv("tar", {"-czf", tarball, "-C", "dist", dist_name})

        print("\n========================================")
        print("Packaged: " .. tarball)
        print("========================================\n")
        print("Distribution structure:")
        print("  " .. dist_name .. "/")
        print("  ├── bin/")
        print("  │   ├── MetaCSSTmain")
        print("  │   ├── MetaCSSTsub")
        print("  │   └── call_vr.py")
        print("  ├── config/")
        print("  │   ├── config.json")
        print("  │   ├── config.toml")
        print("  │   ├── config.yaml")
        print("  │   └── align/          <-- resolves config->align chain")
        print("  │       └── *.align")
        print("  ├── example/")
        print("  │   ├── hv29.fa")
        print("  │   └── out-DGR.gtf")
        print("  ├── requirements.txt")
        print("  └── README")
        print("")
        print("Usage after extraction:")
        print("  ./" .. dist_name .. "/bin/MetaCSSTmain \\")
        print("      -build ./" .. dist_name .. "/config/config.json \\")
        print("      -in     ./" .. dist_name .. "/example/hv29.fa \\")
        print("      -out    result \\")
        print("      -thread 4")
        print("")
        print("  python3 ./" .. dist_name .. "/bin/call_vr.py \\")
        print("      result/raw.gtf \\")
        print("      ./" .. dist_name .. "/example/hv29.fa \\")
        print("      result/final.gtf")
        print("========================================")
    end)
