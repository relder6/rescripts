import argparse

# Parse command-line argument
parser = argparse.ArgumentParser(description="Format ADC timing windows into table")
parser.add_argument("basename", help="Base name for the output file")
args = parser.parse_args()

# Construct output file name
output_filename = f"{args.basename}_timing_window.txt"

input_data = """
hhodo_1x_good_adctdc_diff_time_vs_pmt_pos
hhodo_1x_PosAdcTimeWindowMin= 14.375, 14.375, 14.375, 14.375, 14.375, 14.375, 14.375, 14.375, 14.375, 14.375, 14.375, 14.375, 13.125, 14.375, 14.375, 15.625
hhodo_1x_PosAdcTimeWindowMax= 54.375, 54.375, 54.375, 54.375, 54.375, 54.375, 54.375, 54.375, 54.375, 54.375, 54.375, 54.375, 53.125, 54.375, 54.375, 55.625
hhodo_1x_good_adctdc_diff_time_vs_pmt_neg
hhodo_1x_NegAdcTimeWindowMin= 11.875, 13.125, 13.125, 13.125, 11.875, 13.125, 13.125, 13.125, 13.125, 13.125, 13.125, 13.125, 13.125, 13.125, 13.125, 13.125
hhodo_1x_NegAdcTimeWindowMax= 51.875, 53.125, 53.125, 53.125, 51.875, 53.125, 53.125, 53.125, 53.125, 53.125, 53.125, 53.125, 53.125, 53.125, 53.125, 53.125
hhodo_1y_good_adctdc_diff_time_vs_pmt_pos
hhodo_1y_PosAdcTimeWindowMin= 14.375, 13.125, 14.375, 13.125, 13.125, 13.125, 14.375, 14.375, 14.375, 13.125
hhodo_1y_PosAdcTimeWindowMax= 54.375, 53.125, 54.375, 53.125, 53.125, 53.125, 54.375, 54.375, 54.375, 53.125
hhodo_1y_good_adctdc_diff_time_vs_pmt_neg
hhodo_1y_NegAdcTimeWindowMin= 11.875, 11.875, 13.125, 11.875, 11.875, 11.875, -30.625, 13.125, 13.125, 11.875
hhodo_1y_NegAdcTimeWindowMax= 51.875, 51.875, 53.125, 51.875, 51.875, 51.875, 9.375, 53.125, 53.125, 51.875
hhodo_2x_good_adctdc_diff_time_vs_pmt_pos
hhodo_2x_PosAdcTimeWindowMin= 11.875, 11.875, 13.125, 13.125, 11.875, 11.875, 11.875, 11.875, 11.875, 11.875, 13.125, 11.875, 13.125, 13.125, 13.125, 13.125
hhodo_2x_PosAdcTimeWindowMax= 51.875, 51.875, 53.125, 53.125, 51.875, 51.875, 51.875, 51.875, 51.875, 51.875, 53.125, 51.875, 53.125, 53.125, 53.125, 53.125
hhodo_2x_good_adctdc_diff_time_vs_pmt_neg
hhodo_2x_NegAdcTimeWindowMin= 11.875, 13.125, 13.125, 13.125, 13.125, 13.125, 13.125, 13.125, 13.125, 13.125, 13.125, 13.125, 14.375, 15.625, 13.125, 14.375
hhodo_2x_NegAdcTimeWindowMax= 51.875, 53.125, 53.125, 53.125, 53.125, 53.125, 53.125, 53.125, 53.125, 53.125, 53.125, 53.125, 54.375, 55.625, 53.125, 54.375
hhodo_2y_good_adctdc_diff_time_vs_pmt_pos
hhodo_2y_PosAdcTimeWindowMin= 11.875, 11.875, 13.125, 13.125, 11.875, 11.875, 11.875, 11.875, 11.875, 11.875
hhodo_2y_PosAdcTimeWindowMax= 51.875, 51.875, 53.125, 53.125, 51.875, 51.875, 51.875, 51.875, 51.875, 51.875
hhodo_2y_good_adctdc_diff_time_vs_pmt_neg
hhodo_2y_NegAdcTimeWindowMin= 10.625, 11.875, 13.125, 13.125, 11.875, 11.875, 13.125, 13.125, 11.875, 13.125
hhodo_2y_NegAdcTimeWindowMax= 50.625, 51.875, 53.125, 53.125, 51.875, 51.875, 53.125, 53.125, 51.875, 53.125

"""

# Define which fields we want
window_types = ["PosAdcTimeWindowMin", "PosAdcTimeWindowMax",
                "NegAdcTimeWindowMin", "NegAdcTimeWindowMax"]

# Define which planes to look for
planes = ["1x", "1y", "2x", "2y"]

# Initialize a dictionary to store parsed values
data = {}
for wtype in window_types:
    data[wtype] = {plane: [] for plane in planes}

# Parse each line
for line in input_data.strip().splitlines():
    if '=' not in line:
        continue
    key, values = line.split('=')
    key = key.strip()
    values = [float(v.strip()) for v in values.strip().split(',')]

    parts = key.split('_')
    plane = parts[1]
    wtype = "_".join(parts[2:])

    if wtype in data:
        data[wtype][plane] = values

# Prepare output lines
output_lines = []

for wtype in window_types:
    header = f";{'':<30}{planes[0]:>10}{planes[1]:>10}{planes[2]:>10}{planes[3]:>10}"
    output_lines.append(header)
    output_lines.append(f"hhodo_{wtype:<25} =")

    nrows = max(len(data[wtype][p]) for p in planes)

    for i in range(nrows):
        row = []
        for p in planes:
            try:
                val = data[wtype][p][i]
                row.append(f"{val:>10.2f}")
            except IndexError:
                row.append(f"{0.0:>10.2f}")
        output_lines.append(" " * 32 + "".join(row))
    output_lines.append("")  # Add a blank line between blocks

# Output to file
with open(output_filename, "w") as f:
    for wtype in window_types:
        # Header row
        f.write(f";{'':<30}{planes[0]:>18}{planes[1]:>18}{planes[2]:>18}{planes[3]:>18} \n")
        f.write(f"hhodo_{wtype:<25} = ")

        # Max number of rows for any plane
        nrows = max(len(data[wtype][p]) for p in planes)

        # Data rows
        for i in range(nrows):
            row = []
            for j, p in enumerate(planes):
                try:
                    val = data[wtype][p][i]
                except IndexError:
                    val = 0.0

                # Append formatted value with comma (omit comma on last item)
                comma = "," if j < len(planes) - 1 else ""
                row.append(f"{val:>17.2f}{comma}")
            
            # First line already started with '='
            if i == 0:
                f.write("".join(row) + "\n")
            else:
                f.write(" " * 34 + "".join(row) + "\n")
        f.write("\n")

print(f"✅ Output written to {output_filename}")
