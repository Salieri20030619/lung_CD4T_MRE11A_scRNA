$ErrorActionPreference = 'Stop'
$documentPath = Join-Path $PSScriptRoot 'closeout_delivery\CD4_IKZF2_STAT1_stage_delivery_zh.docx'
$pdfPath = Join-Path $PSScriptRoot 'closeout_delivery\CD4_IKZF2_STAT1_stage_delivery_zh.pdf'
$wordApp = New-Object -ComObject Word.Application
$wordApp.Visible = $false
$wordApp.DisplayAlerts = 0
try {
  $reportDocument = $wordApp.Documents.Open($documentPath, $false, $true)
  $reportDocument.Repaginate()
  $reportDocument.ExportAsFixedFormat($pdfPath, 17)
  $reportDocument.Close(0)
  Write-Output $pdfPath
} finally {
  $wordApp.Quit()
  [System.Runtime.Interopservices.Marshal]::ReleaseComObject($wordApp) | Out-Null
}
