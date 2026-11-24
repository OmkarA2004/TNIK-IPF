import pandas as pd
from openpyxl.styles import Font, PatternFill, Alignment, Border, Side
from openpyxl.utils import get_column_letter
from openpyxl import load_workbook

source = "09_results/tables/TNIK_final_structural_results_table.csv"
output = "09_results/tables/TNIK_final_structural_results_table.xlsx"

df = pd.read_csv(source)
df.to_excel(output, index=False, sheet_name="Structural Results")

wb = load_workbook(output)
ws = wb["Structural Results"]

# Format header
for cell in ws[1]:
    cell.font = Font(bold=True, color="FFFFFF")
    cell.fill = PatternFill("solid", fgColor="1F4E78")
    cell.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)

# Format data cells
for row in ws.iter_rows(min_row=2):
    for cell in row:
        cell.alignment = Alignment(horizontal="center", vertical="center")
        cell.border = Border(bottom=Side(style="thin", color="D9E2F3"))

# Highlight tepotinib row
for cell in ws[3]:
    cell.fill = PatternFill("solid", fgColor="E2F0D9")

# Set readable column widths
for column in ws.columns:
    max_length = max(len(str(cell.value or "")) for cell in column)
    ws.column_dimensions[get_column_letter(column[0].column)].width = min(max_length + 4, 38)

ws.row_dimensions[1].height = 42
ws.freeze_panes = "A2"
ws.auto_filter.ref = ws.dimensions

wb.save(output)
print("Excel workbook created successfully:")
print(output)
