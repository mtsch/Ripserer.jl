using Aqua
using Ripserer

Aqua.test_ambiguities([Ripserer, Base]; exclude=[convert, all, in])
Aqua.test_unbound_args(Ripserer)
Aqua.test_undefined_exports(Ripserer)
Aqua.test_project_extras(Ripserer)
Aqua.test_stale_deps(Ripserer)
Aqua.test_deps_compat(Ripserer)
Aqua.test_project_toml_formatting(Ripserer)
