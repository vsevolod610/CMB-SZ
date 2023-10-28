SCRIPT_PATH=$(dirname "$(readlink -f "$0")")
ROOT="$SCRIPT_PATH/.."
OUTPUT="$ROOT/src/main/output"
TARGET="pic_chain pic_params pic_fit res result*"

cd $OUTPUT
#tar cvzf $ROOT/archive.tar.gz $TARGET
tar cvf $ROOT/archive.tar.gz $TARGET
cd -

# to extract files:
# tar xvf archiv.tar.gz
