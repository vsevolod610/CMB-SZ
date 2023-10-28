SCRIPT_PATH=$(dirname "$(readlink -f "$0")")
SRC="$SCRIPT_PATH/../src"
NULL="/dev/null"

rm $SRC/main/input/szdata/* 2> $NULL && echo "clear main/input/szdata  is DONE!"
rm $SRC/main/input/priors/* 2> $NULL && echo "clear main/input/priors  is DONE!"
rm $SRC/tests/pics/gen/*    2> $NULL && echo "clear tests/pics/gen     is DONE!"
rm $SRC/main/output/pic*/*  2> $NULL && echo "clear main/output/pic*   is DONE!"
rm $SRC/main/output/result* 2> $NULL && echo "clear main/output/result is DONE!"
rm $SRC/main/output/res/*   2> $NULL && echo "clear main/output/res    is DONE!"

#echo -n "folder:"; rm folder/* && echo "is Done!" || echo "is Empty"
