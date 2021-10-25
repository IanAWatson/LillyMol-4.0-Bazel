
cc_binary(
    name = "smu",
    srcs = [
        "smu_main.cc",
        "smu.cc",
        "smu.h",
        "smu_results.cc",
        "smu_results.h",
    ],
    deps = [
         "//Foundational/cmdline:iwcmdline",
         "//Foundational/iwaray:iwaray",
         "//Foundational/iwmisc:iwmisc",
         "//Molecule_Lib:iwmolecule",
         "//Molecule_Lib:moleculeio",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_binary(
    name = "smu_duplicates",
    srcs = [
        "smu_duplicates.cc",
    ],
    deps = [
         "//Foundational/cmdline:iwcmdline",
         "//Foundational/iwaray:iwaray",
         "//Foundational/iwmisc:iwmisc",
         "//Molecule_Lib:iwmolecule",
         "//Molecule_Lib:moleculeio",
    ],
    linkopts = [
        "-lz",
    ],
)


cc_binary(
    name = "smu_rotatable_bond",
    srcs = [
        "smu_rotatable_bond.cc",
    ],
    deps = [
        "//Foundational/accumulator:accumulator",
         "//Foundational/cmdline:iwcmdline",
         "//Foundational/iwaray:iwaray",
         "//Foundational/iwmisc:iwmisc",
         "//Molecule_Lib:iwmolecule",
         "//Molecule_Lib:moleculeio",
    ],
    linkopts = [
        "-lz",
    ],
)
