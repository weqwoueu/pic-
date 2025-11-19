SUBROUTINE My_Label(Code_Name, Code_Description, Code_Developer, Code_Assistant)

IMPLICIT NONE

CHARACTER(*)	Code_Name, Code_Description, Code_Developer, Code_Assistant

CHARACTER*80	TodayDate, OutputDate, StarLine, BlankLine, LabelLine, Label


StarLine	= REPEAT('*',76)
BlankLine	= '*'//REPEAT(' ',74)//'*'

CALL DATE_AND_TIME ( TodayDate )		! date = CCYYMMDD
OutputDate = TodayDate(5:6)//'/'//TodayDate(7:8)//'/'//TodayDate(1:4)

WRITE(6,*) TRIM(StarLine)
WRITE(6,*) TRIM(BlankLine)
CALL Put_Label(TRIM(Code_Name), BlankLine, LabelLine)
WRITE(6,*) TRIM(LabelLine)
CALL Put_Label(TRIM(Code_Description), BlankLine, LabelLine)
WRITE(6,*) TRIM(LabelLine)
CALL Put_Label('Developed by: '//TRIM(Code_Developer), BlankLine, LabelLine)
WRITE(6,*) TRIM(LabelLine)
IF (Code_Assistant/='') THEN
	CALL Put_Label('with the assistance of '//TRIM(Code_Assistant), BlankLine, LabelLine)
	WRITE(6,*) TRIM(LabelLine)
END IF
CALL Put_Label("Computational Advenced Propulsion Lab (CAPLab), Virginia Tech", BlankLine, LabelLine)
WRITE(6,*) TRIM(LabelLine)
CALL Put_Label(TRIM(OutputDate), BlankLine, LabelLine)
WRITE(6,*) TRIM(LabelLine)
WRITE(6,*) TRIM(BlankLine)
WRITE(6,*) TRIM(StarLine)
	
END SUBROUTINE

!	*********************************************************************************************************

SUBROUTINE Put_Label(Label, BlankLine, LabelLine)

CHARACTER*(*)	Label, BlankLine, LabelLine

INTEGER			LabelLen, BlankLineLen, LabelStart, LabelEnd

BlankLineLen	= LEN(BlankLine)
LabelLen		= LEN(TRIM(Label))
LabelStart		= (BlankLineLen - LabelLen)/2 
LabelEnd		= (BlankLineLen + LabelLen)/2 
LabelLine		= BlankLine(1:LabelStart-1)//TRIM(Label)//BlankLine(LabelEnd:BlankLineLen)

END SUBROUTINE

