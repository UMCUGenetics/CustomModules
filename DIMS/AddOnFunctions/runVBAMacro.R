runVBAMacro <- function(dir, dir2, vb_script) {
  
#   dir="E:\\Metabolomics\\Lynne_BSP-2015-08-05\\results\\xls\\"
#   dir2="E:\\Metabolomics\\Lynne_BSP-2015-08-05\\results\\xls_fixed\\"  
  
dir.create(dir2)
files=list.files(dir)

script_1 = paste("Option Explicit
On Error Resume Next
RunExcelMacro

Sub RunExcelMacro() 

Dim xlApp 
Dim xlBook 
Dim xlBook_persenal

Set xlApp = CreateObject(\"Excel.Application\")
xlApp.DisplayAlerts = FALSE
Set xlBook_persenal = xlApp.Workbooks.Open(\"C:\\Users\\awillem8\\AppData\\Roaming\\Microsoft\\Excel\\XLSTART\\PERSONAL.XLSB\", 0, True)
Set xlBook = xlApp.Workbooks.Open(\"", dir, sep="") 

script_2 = "\", 0, True)
xlApp.Run \"PERSONAL.XLSB!FixAndResize\"
xlBook.SaveAs \""

script_3 = "\"
xlBook.Close False 
xlApp.Quit
Set xlBook = Nothing
Set xlBook_persenal = Nothing 
Set xlApp = Nothing 
End Sub"

for (i in 1:length(files)){
  script = paste(script_1, files[i], script_2, dir2, files[i], script_3, sep="")  # dir2, files[i], 
  
  message(script)
  
  fileConn = file("./src/run.vbs")
  writeLines(script, fileConn)
  close(fileConn)
  
  system(vb_script)
}
}

