SCRIPT_PATH=$(dirname "$(readlink -f "$0")")
ROOT="$SCRIPT_PATH/.."

cd $ROOT
python -m venv .venv
source .venv/bin/activate
pip install --upgrade pip
pip install -r requirements.txt
deactivate
cd -
