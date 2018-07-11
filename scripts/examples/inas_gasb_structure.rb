# Geometry

$size = [6, 6]
$latticeConstant = 6.0583
$hsdimension = 3

# Layers

$layers = [ { 
                :name => "InAs Bottom",
                :cation => "In",
                :anion => "As",
                :thickness => 2},
           {
                :name => "GaSb Barrier",
                :cation => "Ga",
                :anion => "Sb",
                :thickness => 2},
           {
                :name => "InAs Top",
                :cation => "In",
                :anion => "As",
                :thickness => 2}]

# Interfaces

$interfaces = [ {
                :layer => 12,
                :material => "InAs",
                :cation => "In", # superfluous since name is given
                :anion => "As", # superfluous since name is given
                :atoms => 128}
            ]

