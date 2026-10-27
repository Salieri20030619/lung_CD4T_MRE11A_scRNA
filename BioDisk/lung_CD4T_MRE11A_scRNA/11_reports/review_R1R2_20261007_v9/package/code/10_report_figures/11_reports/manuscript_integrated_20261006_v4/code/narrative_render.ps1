$ErrorActionPreference = 'Stop'
$documentPath = Join-Path $PSScriptRoot 'narrative_delivery\CD4_IKZF2_STAT1_integrated_manuscript_zh.docx'
$pdfPath = Join-Path $PSScriptRoot 'narrative_delivery\CD4_IKZF2_STAT1_integrated_manuscript_zh.pdf'
$wordApp = New-Object -ComObject Word.Application
$wordApp.Visible = $false
$wordApp.DisplayAlerts = 0
try {
 $reportDocument = $wordApp.Documents.Open($documentPath,$false,$true)
 $reportDocument.Repaginate()
 $reportDocument.ExportAsFixedFormat($pdfPath,17)
 $reportDocument.Close(0)
 Write-Output $pdfPath
} finally {
 $wordApp.Quit()
 [System.Runtime.InteropServices.Marshal]::ReleaseComObject($wordApp) | Out-Null
}
