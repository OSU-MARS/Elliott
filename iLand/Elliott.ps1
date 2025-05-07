$unitTestPath = ([System.IO.Path]::Combine($env:USERPROFILE, "PhD\\iLand\\UnitTests"))
$buildDirectory = ([System.IO.Path]::Combine($unitTestPath, "bin\\x64\\Debug\\net9.0"))
$buildDirectory = ([System.IO.Path]::Combine($unitTestPath, "bin\\x64\\Release\\net9.0"))
Import-Module -Name ([System.IO.Path]::Combine($buildDirectory, "iLand.dll"));
$env:PATH = $env:PATH + (';' + $buildDirectory + '\runtimes\win-x64\native') # for GDAL if GeoTIFF logging is enabled

$elliottPath = ([System.IO.Path]::Combine($env:USERPROFILE, "PhD\\Elliott\\iLand"))

## model instantiation, simulation, and transfer of in memory trajectories to disk
$trajectoryLength = 78 # years
$elliott = Get-Trajectory -Project ([System.IO.Path]::Combine($elliottPath, "Elliott.xml")) -Years $trajectoryLength -Verbose # 2021-2100 -> 78 timesteps past initialization
Write-Trajectory -Trajectory $elliott -ResourceUnitFile ([System.IO.Path]::Combine($elliottPath, "output\\Elliott $trajectoryLength year resource unit trajectories.feather")) -StandFile ([System.IO.Path]::Combine($elliottPath, "output\\Elliott $trajectoryLength year stand trajectories.feather")) -ThreePGFile ([System.IO.Path]::Combine($elliottPath, "output\\Elliott $trajectoryLength year 3-PG.feather")) -Verbose
$elliott.PerformanceCounters
#Write-Trajectory -Trajectory $elliott -IndividualTreeFile ([System.IO.Path]::Combine($elliottPath, "output\\Elliott individual tree trajectories.feather")) -ResourceUnitFile ([System.IO.Path]::Combine($elliottPath, "output\\Elliott resource unit trajectories.feather")) -StandFile ([System.IO.Path]::Combine($elliottPath, "..\\elliottPath\\Elliott stand trajectories.feather")) -ThreePGFile ([System.IO.Path]::Combine($elliottPath, "output\\Elliott 3-PG.feather"))
