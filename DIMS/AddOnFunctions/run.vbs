Option Explicit
On Error Resume Next
RunExcelMacro

Sub RunExcelMacro() 

Dim xlApp 
Dim xlBook 
Dim xlBook_persenal

Set xlApp = CreateObject("Excel.Application")
xlApp.DisplayAlerts = FALSE
Set xlBook_persenal = xlApp.Workbooks.Open("C:\Users\awillem8\AppData\Roaming\Microsoft\Excel\XLSTART\PERSONAL.XLSB", 0, True)
Set xlBook = xlApp.Workbooks.Open("C:\Users\mkerkho7\testmap\xls\P181_Aberant_HMDB.xlsx", 0, True)
xlApp.Run "PERSONAL.XLSB!FixAndResize"
xlBook.SaveAs "C:\Users\mkerkho7\testmap\results\P181_Aberant_HMDB.xlsx"
xlBook.Close False 
xlApp.Quit
Set xlBook = Nothing
Set xlBook_persenal = Nothing 
Set xlApp = Nothing 
End Sub
