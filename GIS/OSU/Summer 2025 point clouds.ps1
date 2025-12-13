$buildDirectory = ([System.IO.Path]::Combine($env:USERPROFILE, "PhD\tools\Clouds\UnitTests\bin\Debug\net9.0"))
#$buildDirectory = ([System.IO.Path]::Combine($env:USERPROFILE, "PhD\tools\Clouds\UnitTests\bin\Release\net9.0"))
$env:PATH = $env:PATH + (';' + $buildDirectory + '\runtimes\win-x64\native') # for GDAL

Import-Module -Name ([System.IO.Path]::Combine($buildDirectory, "Clouds.dll"))

$projectDirectory = "D:\Elliott\GIS\OSU" # update as needed

# s03990w06810
# 207797.650 m, 121855.642 m, 474.421 m
# 207789.110 m, 121851.372 m, 474.691 m, 18°
Register-Clouds -Las "$projectDirectory\0717_38 classified.las" -Lat 43.56772995 -Long -123.94454956 -Z (3.28084 * 474.691) -NudgeY (-3.28084 * 8.54) -NudgeX (-3.28084 * 4.27) -RotationXY 18 -HorizontalEpsg 6557 -VerticalEpsg 8228 -FallbackDate "2025-07-17" -RepairBounds -RepairReturn