# list source files
set (FILES
  Evaluator
  InterfaceGenerator
  LatticeGenerator
  Optimizer
  )

if (VTK_FOUND)
  set (FILES
    ${FILES}
    Extractor
    Plot
    Visualize
    )
endif()
