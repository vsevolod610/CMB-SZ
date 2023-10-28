SCRIPT_PATH=$(dirname "$(readlink -f "$0")")
ROOT="$SCRIPT_PATH/.."
SRC="$ROOT/src"
NULL="/dev/null"


cd $SRC/main
#python check.py
python config.py
read -p "Continue? (Y/N): " confirm && [[ $confirm == [yY] || $confirm == [yY][eE][sS] ]] || exit 1
echo "Start: $(date)" >> $ROOT/logs/log.txt
rm $SRC/main/output/result* 2> $NULL && echo "clear main/output/result is DONE!"
python create_data.py || exit 1
python start_analyze.py || exit 1
python estimate.py > $SRC/main/output/res/T0.txt || exit 1
cd -
echo "Finish: $(date)" >> $ROOT/logs/log.txt
