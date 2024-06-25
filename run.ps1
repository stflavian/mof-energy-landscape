
$SOFTWARE_PATH = "full/path/to/this/folder$

$SOFTWARE = "$SOFTWARE_PATH\src\main.jl"
$ENVIRONMENT = "$SOFTWARE_PATH\."

julia --project=$ENVIRONMENT $SOFTWARE
