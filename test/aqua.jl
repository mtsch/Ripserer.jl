using Aqua
using Ripserer

if VERSION < v"1.8-DEV"
    # Note: will have to change Project.toml for 1.8
    Aqua.test_all(Ripserer; ambiguities=false)
end
