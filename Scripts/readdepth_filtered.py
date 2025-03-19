import os


folder_path = './'
depth_files = [f for f in os.listdir(folder_path) if f.endswith("_100readdepth.txt")]

for depth_file in depth_files:

    file_path = os.path.join(folder_path, depth_file)

    filtered_lines = []
    with open(file_path, "r", encoding='utf-8') as file:
        for line in file:
            parts = line.strip().split("\t")
            try:
                if len(parts) >= 4 and float(parts[3]) == 0:
                    filtered_lines.append("\t".join(parts[:4]))
            except ValueError:
                continue

    output_file_name = depth_file.replace("_100readdepth.txt", "_filtered100readdepth.txt")
    output_file_path = os.path.join(folder_path, output_file_name)
    with open(output_file_path, "w", encoding='utf-8') as output_file:
        for line in filtered_lines:
            output_file.write(line + "\n")

print("Filtering completed.")
