SCRIPT_PATH=$(dirname "$(readlink -f "$0")")
SRC="$SCRIPT_PATH/../src"

cd "$SRC/tests"
SHOW=$(ls | grep "^show.*.py")

for test in $SHOW;
do
    echo "Run $test ..."
    python $test | tee res/${test%.py}.txt || exit 1
    echo " "
done

cd -
