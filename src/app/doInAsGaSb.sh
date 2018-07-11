
./LatticeGenerator --config ../xml/latticeConfig.xml
./LatticeGenerator --config ../xml/latticeConfig.xml --material-name GaSb -i ../data/InAs_GaSb/lattice.xml --xyz ../data/InAs_GaSb/step2.xyz
./LatticeGenerator --config ../xml/latticeConfig.xml -i ../data/InAs_GaSb/lattice.xml --xyz ../data/InAs_GaSb/step3.xyz
./InterfaceGenerator --config ../xml/interfaceConfig.xml
#./Optimizer --config ../xml/optimizerConf.xml
