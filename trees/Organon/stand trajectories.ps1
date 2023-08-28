Set-Location -Path ([System.IO.Path]::Combine($env:USERPROFILE, "PhD\\Elliott\\trees\\Organon"))
$buildDirectory = ([System.IO.Path]::Combine($env:USERPROFILE, "PhD\\SEEM\\UnitTests\\bin\\x64\\Debug\\net7.0-windows"))
#$buildDirectory = ([System.IO.Path]::Combine($env:USERPROFILE, "PhD\\SEEM\\UnitTests\\bin\\x64\\Release\\net7.0-windows"))
#$buildDirectory = ([System.IO.Path]::Combine($env:USERPROFILE, "PhD\\SEEM\\Seem\\bin\\win-x64\\Release\\net7.0\\win-x64\publish"))
Import-Module -Name ([System.IO.Path]::Combine($buildDirectory, "Seem.dll"))

$financial = Get-FinancialScenarios -Xlsx ([System.IO.Path]::Combine((Get-Location), "financial scenarios.xlsx")) -XlsxSheet "parameterization"

# 142 cruise stands assigned for intensive management
$intensiveStands = Get-CruisedStands -Model OrganonSWO -TreesSheet "intensiveTrees" -Xlsx ([System.IO.Path]::Combine((Get-Location), "Elliott Organon cruise records 2015-16.xlsx"))
$allTrajectories = New-Object "System.Collections.Generic.List[Mars.Seem.Heuristics.HeuristicStandTrajectories[Mars.Seem.Heuristics.PrescriptionParameters]]"

for ($standIndex = 0; $standIndex -lt $intensiveStands.Stands.Count; ++$standIndex)
{
    Write-Host "Stand $($intensiveStands.Stands[$standIndex].Name)..."
    $standTrajectories = Optimize-Prescription -Stand $intensiveStands.Stands[$standIndex] -TreeModel $intensiveStands.OrganonVariant.TreeModel -Enumerate -FirstThinAge (-1, 35, 40, 45, 50, 55) -RotationAge (50, 55, 60, 65, 70, 75) -FromAbovePercentageUpperLimit 0 -FromBelowPercentageUpperLimit 0 -MinimumIntensity 20 -MaximumIntensity 40 -DefaultStep 5 -MinimumStep 5
    $allTrajectories.Add($standTrajectories)
}

Write-SilviculturalTrajectories -Trajectories $allTrajectories -StartYear 2016 -FilePath ([System.IO.Path]::Combine((Get-Location), "Elliott intensive prescriptions.csv")) -NoCarbon
Write-SilviculturalTrajectories -Trajectories $allTrajectories -StartYear 2016 -FilePath ([System.IO.Path]::Combine((Get-Location), "Elliott intensive prescriptions.feather")) -NoCarbon

# no management baseline: site indices from GIS
$stands = Get-CruisedStands -Model OrganonSWO -Xlsx ([System.IO.Path]::Combine((Get-Location), "Elliott Organon cruise records 2015-16.xlsx"))
$standTrajectories = Get-StandTrajectories -Stands $stands
Write-StandTrajectories -Trajectories $standTrajectories -Financial $financial -StartYear 2016 -FilePath ([System.IO.Path]::Combine((Get-Location), "Elliott stand trajectories 2016-2116.feather")) -NoCarbon
Write-StandTrajectories -Trajectories $standTrajectories -Financial $financial -StartYear 2016 -FilePath ([System.IO.Path]::Combine((Get-Location), "Elliott stand trajectories 2016-2116.csv")) -NoCarbon

# no management baseline: site index sensitvity
# Using % in a file name crashes R, so write files with "percent" instead
$stands90 = Get-CruisedStands -Model OrganonSWO -Xlsx ([System.IO.Path]::Combine((Get-Location), "Elliott Organon cruise records 90% site index.xlsx"))
$standTrajectories90 = Get-StandTrajectories -Stands $stands90
Write-StandTrajectories -Trajectories $standTrajectories90 -Financial $financial -CsvFile ([System.IO.Path]::Combine((Get-Location), "Elliott stand trajectories 90 percent site index.csv")) -NoCarbon

$stands80 = Get-CruisedStands -Model OrganonSWO -Xlsx ([System.IO.Path]::Combine((Get-Location), "Elliott Organon cruise records 80% site index.xlsx"))
$standTrajectories80 = Get-StandTrajectories -Stands $stands80
Write-StandTrajectories -Trajectories $standTrajectories80 -Financial $financial -CsvFile ([System.IO.Path]::Combine((Get-Location), "Elliott stand trajectories 80 percent site index.csv")) -NoCarbon