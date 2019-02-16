## Sequence

```mermaid
sequenceDiagram
  participant user as User
  participant atoms as BondedAtoms
  participant writer as LammpsWriter
  user ->> atoms : create
  atoms -->> user : 'atoms'
  note over user, atoms : setting bonds and types, ..., etc.
  user ->> writer : create with 'atoms'
  user ->> writer : call set_atom_data( )
  user ->> writer : call set_bond_types( )
  user ->> writer : call set_angle_types( )
  user ->> writer : call set_dihedral_types( )
  user ->> writer : call set_improper_types( )
  user ->> writer : call write_lammps_data( )
  note over writer : write a data file!
```

## Class

```puml
package ase {
  class Atoms
  package calculators.lammpsrun {
    class Prism
  }
}

class BondedAtoms {

}

note left: Has ASE's functionalities

class ExtendedPrism {

}

class LammpsWriter {

}

class LammpsAtoms {

}

class LammpsTopology {

}

class LammpsBonds {
}

class LammpsAngles {
}

class LammpsDihedrals {
}

class LammpsImpropers {

}

class LammpsSpecialBonds {

}

Atoms <|-- BondedAtoms
Prism <|-- ExtendedPrism

LammpsTopology <|-- LammpsBonds
LammpsTopology <|-- LammpsAngles
LammpsTopology <|-- LammpsDihedrals
LammpsTopology <|-- LammpsImpropers

LammpsWriter "1" o-- "1" LammpsAtoms
LammpsWriter "1" o-- "4" LammpsTopology
LammpsWriter "1" o-- "1" LammpsSpecialBonds

BondedAtoms <.. LammpsWriter
ExtendedPrism <.. LammpsWriter
```
