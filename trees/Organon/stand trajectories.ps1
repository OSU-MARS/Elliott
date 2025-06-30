Set-Location -Path ([System.IO.Path]::Combine($env:USERPROFILE, "PhD\\Elliott\\trees\\Organon"))
$buildDirectory = ([System.IO.Path]::Combine($env:USERPROFILE, "PhD\\SEEM\\UnitTests\\bin\\x64\\Debug\\net9.0-windows10.0.19041.0"))
$buildDirectory = ([System.IO.Path]::Combine($env:USERPROFILE, "PhD\\SEEM\\UnitTests\\bin\\x64\\Release\\net9.0-windows10.0.19041.0"))
#$buildDirectory = ([System.IO.Path]::Combine($env:USERPROFILE, "PhD\\SEEM\\Seem\\bin\\win-x64\\Release\\net9.0-windows10.0.19041.0\\publish"))
Import-Module -Name ([System.IO.Path]::Combine($buildDirectory, "Seem.dll"))

$financial = Get-FinancialScenarios -Xlsx ([System.IO.Path]::Combine((Get-Location), "financial scenarios.xlsx")) -XlsxSheet "parameterization"


# no management baseline: site indices from GIS
$stands = Get-CruisedStands -Model OrganonSWO -Xlsx ([System.IO.Path]::Combine((Get-Location), "Elliott Organon cruise records 2015-16 v2.xlsx"))
$simulationTime = Measure-Command { $standTrajectories = Get-StandTrajectories -Stands $stands -Years 100 }

Write-StandTrajectories -Trajectories $standTrajectories -Financial $financial -StartYear 2016 -FilePath ([System.IO.Path]::Combine((Get-Location), "Elliott stand trajectories 2016-2116.feather")) -NoCarbon
#Write-StandTrajectories -Trajectories $standTrajectories -Financial $financial -StartYear 2016 -FilePath ([System.IO.Path]::Combine((Get-Location), "Elliott stand trajectories 2016-2116.csv")) -NoCarbon

Write-TreeList -Trajectories $standTrajectories -StartYear 2016 -FilePath ([System.IO.Path]::Combine((Get-Location), "Elliott tree lists 2016-2116.feather"))
#Write-TreeList -Trajectories $standTrajectories -StartYear 2016 -FilePath ([System.IO.Path]::Combine((Get-Location), "Elliott tree lists 2016-2116.csv"))


# cruised plantations: LEV maximization
$plantations = Get-CruisedStands -Model OrganonSWO -TreesSheet "plantationTrees" -Xlsx ([System.IO.Path]::Combine((Get-Location), "Elliott Organon cruise records 2015-16 v2.xlsx"))
$thinAges = (-1, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85)
$rotationAges = (35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 115, 120, 125)

$allTrajectories = New-Object "System.Collections.Generic.List[Mars.Seem.Heuristics.HeuristicStandTrajectories[Mars.Seem.Heuristics.PrescriptionParameters]]"
# 44 of 418 plantation stands trigger modeling limitations, 41 of which are marginal for thinning and thus of little interest
# no thin found, hits null as unthinned prescription isn't set: 34, 84, 176, 227, 228, 351, 515, 546, 552, 604, 993, 1007, 1436, 1664, 1667, 1754, 1967, 2011, 2160, 2161, 2359, 2389, 2404, 2407, 2411, 2414, 2416, 2421, 2432, 2462
# mismatch in volume calculations has forwarder forwarding zero merch m3: 50, 58, 72, 1045, 2373, 2427, 2445, 2461, 2531
# mismatch in volume calculations has zero merch m3 yarded: 52, 2450
# duplicate tags trigger compaction asserts: 2437, 2465
# intermittent coordinate not found: 306
# (Also 546, 2455 generate Organon diameter growth multipliers > 1.5.)
$standsToBypassForNow = ("34", "50", "52", "58", "72", "84", "176", "227", "228", "351", "515", "546", "552", "604", "993", "1007", "1045", "1436", "1664", "1667", "1754", "1967", "2011", "2160", "2161", "2359", "2373", "2389", "2404", "2407", "2411", "2414", "2416", "2421", "2427", "2432", "2437", "2445", "2450", "2461", "2462", "2465", "2531")
for ($standIndex = 0; $standIndex -lt $plantations.Stands.Count; ++$standIndex)
{
	$stand = $plantations.Stands[$standIndex]
	if ($standsToBypassForNow.Contains($stand.Name))

	{
		continue
	}
	
    Write-Host "Stand $($plantations.Stands[$standIndex].Name)..."
    $standTrajectories = Optimize-Prescription -Stand $stand -Financial $financial -TreeModel $plantations.OrganonVariant.TreeModel -Enumerate -FirstThinAge $thinAges -RotationAge $rotationAges -FromAbovePercentageUpperLimit 0 -FromBelowPercentageUpperLimit 50 -MinimumIntensity 10 -MaximumIntensity 50 -DefaultStep 5 -MinimumStep 1.25
    $allTrajectories.Add($standTrajectories)
}

# TODO: .feather size is misestimated at 30 GB instead of 4.1 GB
Write-SilviculturalTrajectories -LimitGB 30 -Trajectories $allTrajectories -StartYear 2016 -FilePath ([System.IO.Path]::Combine((Get-Location), "Elliott plantation prescriptions max LEV.feather")) -NoCarbon
#Write-SilviculturalTrajectories -Trajectories $allTrajectories -StartYear 2016 -FilePath ([System.IO.Path]::Combine((Get-Location), "Elliott plantation prescriptions max LEV.csv")) -NoCarbon

# cruised plantations: NPV maximization
$allTrajectories = New-Object "System.Collections.Generic.List[Mars.Seem.Heuristics.HeuristicStandTrajectories[Mars.Seem.Heuristics.PrescriptionParameters]]"
for ($standIndex = 0; $standIndex -lt $plantations.Stands.Count; ++$standIndex)
{
	$stand = $plantations.Stands[$standIndex]
	if ($standsToBypassForNow.Contains($stand.Name))
	{
		continue
	}
	
    Write-Host "Stand $($plantations.Stands[$standIndex].Name)..."
    $standTrajectories = Optimize-Prescription -Stand $plantations.Stands[$standIndex] -Financial $financial -TimberObjective NetPresentValue -TreeModel $plantations.OrganonVariant.TreeModel -Enumerate -FirstThinAge $thinAges -RotationAge $rotationAges -FromAbovePercentageUpperLimit 0 -FromBelowPercentageUpperLimit 40 -MinimumIntensity 10 -MaximumIntensity 45 -DefaultStep 5 -MinimumStep 2.5
    $allTrajectories.Add($standTrajectories)
}

Write-SilviculturalTrajectories -Trajectories $allTrajectories -StartYear 2016 -FilePath ([System.IO.Path]::Combine((Get-Location), "Elliott intensive prescriptions max NPV.csv")) -NoCarbon
Write-SilviculturalTrajectories -Trajectories $allTrajectories -StartYear 2016 -FilePath ([System.IO.Path]::Combine((Get-Location), "Elliott intensive prescriptions max NPV.feather")) -NoCarbon


# no management baseline: site index sensitvity
# Using % in a file name crashes R, so write files with "percent" instead
$stands90 = Get-CruisedStands -Model OrganonSWO -Xlsx ([System.IO.Path]::Combine((Get-Location), "Elliott Organon cruise records 90% site index.xlsx"))
$standTrajectories90 = Get-StandTrajectories -Stands $stands90
Write-StandTrajectories -Trajectories $standTrajectories90 -Financial $financial -CsvFile ([System.IO.Path]::Combine((Get-Location), "Elliott stand trajectories 90 percent site index.csv")) -NoCarbon

$stands80 = Get-CruisedStands -Model OrganonSWO -Xlsx ([System.IO.Path]::Combine((Get-Location), "Elliott Organon cruise records 80% site index.xlsx"))
$standTrajectories80 = Get-StandTrajectories -Stands $stands80
Write-StandTrajectories -Trajectories $standTrajectories80 -Financial $financial -CsvFile ([System.IO.Path]::Combine((Get-Location), "Elliott stand trajectories 80 percent site index.csv")) -NoCarbon