$ErrorActionPreference = 'Stop'
$documentPath = Join-Path $PSScriptRoot 'graphics_delivery\CD4_IKZF2_STAT1_integrated_manuscript_zh.docx'
$pdfPath = Join-Path $PSScriptRoot 'graphics_delivery\CD4_IKZF2_STAT1_integrated_manuscript_zh.pdf'
$wordApp = New-Object -ComObject Word.Application
$wordApp.Visible = $false
$wordApp.DisplayAlerts = 0
try {
 Write-Output 'OPENING_REPORT'
 $reportDocument = $wordApp.Documents.Open($documentPath,$false,$true)
 Write-Output 'EXPORTING_PDF'
 $reportDocument.ExportAsFixedFormat($pdfPath,17)
 Write-Output 'PDF_EXPORTED'
 $reportDocument.Close(0)
 Write-Output $pdfPath
} finally {
 $wordApp.Quit()
 [System.Runtime.InteropServices.Marshal]::ReleaseComObject($wordApp) | Out-Null
}
