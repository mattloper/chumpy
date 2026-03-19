#!/usr/bin/env bash
set -e

echo "════════════════════════════════════════"
echo "  Testing chumpy build and install"
echo "════════════════════════════════════════"
echo ""

# Clean previous artifacts
echo "🧹 Cleaning previous builds..."
rm -rf dist/ build/ *.egg-info test-pip test-uv test-wheel

echo ""
echo "════════════════════════════════════════"
echo "  Test 1: pip install ."
echo "════════════════════════════════════════"
echo ""

python3 -m venv test-pip
source test-pip/bin/activate
pip install . -q
python -c "import chumpy; print(f'✓ Test 1 PASSED: chumpy {chumpy.__version__} via pip install .')"
deactivate
rm -rf test-pip

echo ""
echo "════════════════════════════════════════"
echo "  Test 2: uv pip install ."
echo "════════════════════════════════════════"
echo ""

uv venv test-uv
source test-uv/bin/activate
uv pip install .
python -c "import chumpy; print(f'✓ Test 2 PASSED: chumpy {chumpy.__version__} via uv pip install .')"
deactivate
rm -rf test-uv

echo ""
echo "════════════════════════════════════════"
echo "  Test 3: python -m build"
echo "════════════════════════════════════════"
echo ""

pip install build -q
python -m build

echo ""
echo "Built artifacts:"
ls -lh dist/

echo ""
echo "Testing wheel installation..."
python3 -m venv test-wheel
source test-wheel/bin/activate

# Install the wheel (use wildcard to match whatever version was built)
WHEEL_FILE=$(ls dist/*.whl | head -n 1)
echo "Installing: $WHEEL_FILE"
pip install "$WHEEL_FILE" -q

python -c "import chumpy; print(f'✓ Test 3 PASSED: chumpy {chumpy.__version__} from wheel')"
deactivate
rm -rf test-wheel

echo ""
echo "════════════════════════════════════════"
echo "  ✅ All tests PASSED!"
echo "════════════════════════════════════════"
echo ""
echo "Ready to:"
echo "  1. Submit PR to original chumpy repo"
echo "  2. Publish to PyPI (optional)"
echo ""
