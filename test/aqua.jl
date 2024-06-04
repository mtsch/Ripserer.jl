using Aqua
using Ripserer

Aqua.test_all(Ripserer; ambiguities=false, piracies=false)
Aqua.test_piracies(Ripserer; treat_as_own=PersistenceDiagram)
