$baseUrl = "http://www.flooxs.ece.ufl.edu"
$startPage = "/index.php/Main_Page"
$outDir = "flooxs_docs"
$maxParallel = 10

New-Item -ItemType Directory -Force -Path $outDir | Out-Null

# Phase 1: Discover all pages (sequential but fast - just parsing)
Write-Host "Discovering pages..."
$visited = @{}
$allPages = [System.Collections.ArrayList]::new()
$queue = [System.Collections.Queue]::new()
$queue.Enqueue($startPage)

while ($queue.Count -gt 0) {
    $page = $queue.Dequeue()
    if ($visited.ContainsKey($page)) { continue }
    $visited[$page] = $true
    [void]$allPages.Add($page)

    $url = $baseUrl + $page
    try {
        $response = Invoke-WebRequest -Uri $url -UseBasicParsing -TimeoutSec 10
        $links = [regex]::Matches($response.Content, 'href="(/index\.php/[^"#]+)"')
        foreach ($link in $links) {
            $href = $link.Groups[1].Value
            if ($href -match 'Special:|action=|title=|oldid=|Talk:|User:|File:|redirect') { continue }
            if (-not $visited.ContainsKey($href)) {
                $queue.Enqueue($href)
            }
        }
        Write-Host "`rFound $($allPages.Count) pages..." -NoNewline
    }
    catch { }
}

Write-Host "`nFound $($allPages.Count) total pages. Downloading in parallel..."

# Phase 2: Download all pages in parallel
$allPages | ForEach-Object -ThrottleLimit $maxParallel -Parallel {
    $page = $_
    $baseUrl = $using:baseUrl
    $outDir = $using:outDir

    $url = $baseUrl + $page
    $filename = $page -replace '/index.php/', '' -replace '/', '_' -replace '\?.*', ''
    if ($filename -eq '') { $filename = 'Main_Page' }
    $filepath = Join-Path $outDir "$filename.html"

    try {
        $response = Invoke-WebRequest -Uri $url -UseBasicParsing -TimeoutSec 30
        $response.Content | Out-File -FilePath $filepath -Encoding UTF8
        Write-Host "OK: $filename"
    }
    catch {
        Write-Host "FAIL: $filename"
    }
}

Write-Host "`nDone!"
