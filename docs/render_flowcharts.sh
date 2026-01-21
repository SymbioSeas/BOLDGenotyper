#!/bin/bash
#
# BOLDGenotyper Flow Chart Renderer
# Generates high-quality SVG and PNG versions of workflow diagrams
#
# Requirements:
#   - Node.js installed
#   - Mermaid CLI: npm install -g @mermaid-js/mermaid-cli
#
# Usage:
#   bash render_flowcharts.sh

set -e  # Exit on error

# Color output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}========================================${NC}"
echo -e "${BLUE}BOLDGenotyper Flow Chart Renderer${NC}"
echo -e "${BLUE}========================================${NC}"
echo ""

# Check if mmdc is installed
if ! command -v mmdc &> /dev/null; then
    echo -e "${RED}Error: Mermaid CLI (mmdc) not found!${NC}"
    echo ""
    echo "Please install Mermaid CLI first:"
    echo "  npm install -g @mermaid-js/mermaid-cli"
    echo ""
    exit 1
fi

echo -e "${GREEN}✓ Mermaid CLI found${NC}"
echo ""

# Get script directory (where this script is located)
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
cd "$SCRIPT_DIR"

echo -e "${BLUE}Working directory: ${SCRIPT_DIR}${NC}"
echo ""

# ========================================
# Render Simplified Flow Chart (Main Text)
# ========================================
echo -e "${YELLOW}[1/4] Rendering simplified flow chart to SVG...${NC}"
mmdc -i flowchart_main_text_simplified.mmd \
     -o flowchart_main_text_simplified.svg \
     -b transparent

if [ -f "flowchart_main_text_simplified.svg" ]; then
    SIZE=$(du -h flowchart_main_text_simplified.svg | cut -f1)
    echo -e "${GREEN}✓ Created: flowchart_main_text_simplified.svg (${SIZE})${NC}"
else
    echo -e "${RED}✗ Failed to create SVG${NC}"
    exit 1
fi

echo -e "${YELLOW}[2/4] Rendering simplified flow chart to PNG (300 DPI)...${NC}"
mmdc -i flowchart_main_text_simplified.mmd \
     -o flowchart_main_text_simplified.png \
     -w 3000 \
     -H 4000 \
     -b white

if [ -f "flowchart_main_text_simplified.png" ]; then
    SIZE=$(du -h flowchart_main_text_simplified.png | cut -f1)
    echo -e "${GREEN}✓ Created: flowchart_main_text_simplified.png (${SIZE})${NC}"
else
    echo -e "${RED}✗ Failed to create PNG${NC}"
    exit 1
fi

echo ""

# ========================================
# Render Detailed Flow Chart (Supplemental)
# ========================================
echo -e "${YELLOW}[3/4] Rendering detailed flow chart to SVG...${NC}"
mmdc -i flowchart_supplemental_detailed.mmd \
     -o flowchart_supplemental_detailed.svg \
     -b transparent

if [ -f "flowchart_supplemental_detailed.svg" ]; then
    SIZE=$(du -h flowchart_supplemental_detailed.svg | cut -f1)
    echo -e "${GREEN}✓ Created: flowchart_supplemental_detailed.svg (${SIZE})${NC}"
else
    echo -e "${RED}✗ Failed to create SVG${NC}"
    exit 1
fi

echo -e "${YELLOW}[4/4] Rendering detailed flow chart to PNG (300 DPI)...${NC}"
mmdc -i flowchart_supplemental_detailed.mmd \
     -o flowchart_supplemental_detailed.png \
     -w 3500 \
     -H 8000 \
     -b white

if [ -f "flowchart_supplemental_detailed.png" ]; then
    SIZE=$(du -h flowchart_supplemental_detailed.png | cut -f1)
    echo -e "${GREEN}✓ Created: flowchart_supplemental_detailed.png (${SIZE})${NC}"
else
    echo -e "${RED}✗ Failed to create PNG${NC}"
    exit 1
fi

echo ""
echo -e "${BLUE}========================================${NC}"
echo -e "${GREEN}✓ All flow charts rendered successfully!${NC}"
echo -e "${BLUE}========================================${NC}"
echo ""
echo "Output files:"
echo "  Main Text Figure (Simplified):"
echo "    - flowchart_main_text_simplified.svg"
echo "    - flowchart_main_text_simplified.png"
echo ""
echo "  Supplemental Figure (Detailed):"
echo "    - flowchart_supplemental_detailed.svg"
echo "    - flowchart_supplemental_detailed.png"
echo ""
echo -e "${YELLOW}Tip:${NC} Use SVG for LaTeX/vector graphics, PNG for Word/PowerPoint"
echo ""
