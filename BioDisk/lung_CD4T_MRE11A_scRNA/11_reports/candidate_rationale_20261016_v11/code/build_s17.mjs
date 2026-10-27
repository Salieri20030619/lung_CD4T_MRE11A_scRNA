import fs from 'node:fs/promises';
import { Workbook, SpreadsheetFile } from '@oai/artifact-tool';
const root='C:/codex_test/tmp/silicosis_v11';
const data=JSON.parse(await fs.readFile(root+'/workbook_data.json','utf8'));
const wb=Workbook.create();
function letter(n){let s='';for(n++;n>0;n=Math.floor((n-1)/26))s=String.fromCharCode(65+(n-1)%26)+s;return s;}
const widths={基因:18,代表状态:12,状态:12,log2FC:12,原始P:14,原BH:14,logCPM:12,候选基因:18,同方向效应更低的候选:36,历史记录:44,文献具体发现:50,研究场景:28,'当前未优先/保留依据':58,当前处置建议:24,再次考虑的条件:44,文献:65,来源链接:65,对象:28,历史明确记录:66,记录编号:14,解释边界:60,原始对话:65,名义基因绝对效应排名:20,同方向且效应更低的候选:40,不分方向效应更低的候选:48,详细文献评估:42,项目:26,说明:110,该状态进入前20展示:24,检验状态:22};
for(let k=0;k<data.sheets.length;k++){
 const d=data.sheets[k],sh=wb.worksheets.add(d.name);sh.showGridLines=false;
 const n=d.rows.length,m=d.columns.length,last=letter(m-1);
 sh.getRange('A2').values=[['S17 '+d.name]];sh.getRange('A2').format.font={name:'Microsoft YaHei',size:14,bold:true,color:'#183448'};
 sh.getRange(`A4:${last}${n+4}`).values=[d.columns,...d.rows];
 const all=sh.getRange(`A4:${last}${n+4}`);all.format.font={name:'Microsoft YaHei',size:10,color:'#243746'};all.format.wrapText=true;all.format.verticalAlignment='top';
 sh.getRange(`A4:${last}4`).format={fill:'#183448',font:{name:'Microsoft YaHei',bold:true,color:'#ffffff',size:10},wrapText:true,rowHeight:42,horizontalAlignment:'center',verticalAlignment:'center'};
 for(let c=0;c<m;c++){
  const range=sh.getRange(`${letter(c)}5:${letter(c)}${n+4}`);range.format.columnWidth=widths[d.columns[c]]||22;
  if(['log2FC','logCPM'].includes(d.columns[c])){range.setNumberFormat('0.0000;[Red]-0.0000');range.format.horizontalAlignment='right';}
  if(['原始P','原BH'].includes(d.columns[c])){range.setNumberFormat('0.0000E+00');range.format.horizontalAlignment='right';}
 }
 sh.getRange(`A5:${last}${n+4}`).format.rowHeight=d.name==='逐基因说明'?72:d.name==='历史筛选依据'?42:d.name==='使用说明'?44:d.name==='完整效应比较'?34:28;
 sh.freezePanes.freezeRows(4);sh.freezePanes.freezeColumns(d.name==='逐基因说明'?2:1);
 const table=sh.tables.add(`A4:${last}${n+4}`,true,'S17Table'+k);table.showFilterButton=true;
 if(d.name==='逐基因说明')sh.tabColor='#183448';
}
wb.recalculate();
const inspection=await wb.inspect({kind:'sheet',include:'id,name',maxChars:1500});console.log(inspection.ndjson||JSON.stringify(inspection));
await fs.mkdir(root+'/previews',{recursive:true});
for(const d of data.sheets){
 const ranges=d.name==='逐基因说明'?['A2:F8','G4:L7']:d.name==='历史筛选依据'?['A2:D9']:d.name==='使用说明'?['A2:B10']:['A2:F11'];
 for(let i=0;i<ranges.length;i++){
  const preview=await wb.render({sheetName:d.name,range:ranges[i],scale:1,format:'png'});
  await fs.writeFile(root+`/previews/${d.name}_${i}.png`,new Uint8Array(await preview.arrayBuffer()));
 }
}
const file=await SpreadsheetFile.exportXlsx(wb);await file.save('C:/codex_test/output/silicosis_v11/supplement/S17_candidate_selection_rationale.xlsx');
console.log('XLSX EXPORTED');
