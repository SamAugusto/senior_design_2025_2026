import tkinter as tk
from tkinter import filedialog, messagebox
from tkinter import ttk
import os
import csv
from PIL import Image, ImageTk, ImageDraw, ImageFont

# add analysis code in the lines 350 later one when finalized - still developing other functionalities of the code
class EnhancedFileApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Automated Detection of Telomeric Structural Variations")
        self.root.geometry("1400x900")
        
        # Variables to store file paths
        self.csv_file_path = tk.StringVar()
        self.txt_file_path = tk.StringVar()
        
        # Variables for table options
        self.option_vars = {}
        self.option_data = {
            "Option 1": [],
            "Option 2": [],
            "Option 3": [],
            "Option 4": []
        }
        
        # Image variable
        self.output_image = None
        
        # Create the GUI components
        self.create_widgets()
    
    def create_widgets(self):
        # Main container
        main_frame = ttk.Frame(self.root, padding="10")
        main_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        
        # Configure grid weights
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(0, weight=1)
        main_frame.columnconfigure(0, weight=1)
        main_frame.columnconfigure(1, weight=1)
        main_frame.rowconfigure(2, weight=1)
        main_frame.rowconfigure(3, weight=1)
        # Title
        title_label = ttk.Label(main_frame, text="Automated Detection of Telomeric Structural Variations", 
                               font=('Arial', 18, 'bold'))
        title_label.grid(row=0, column=0, columnspan=2, pady=(0, 15))
        
        # ===== LEFT COLUMN =====
        left_frame = ttk.Frame(main_frame)
        left_frame.grid(row=1, column=0, rowspan=3, sticky=(tk.W, tk.E, tk.N, tk.S), 
                       padx=(0, 10))
        left_frame.columnconfigure(0, weight=1)
        left_frame.rowconfigure(2, weight=1)
        left_frame.rowconfigure(4, weight=1)
        
        # CSV File Input Section
        csv_section = ttk.LabelFrame(left_frame, text="Result File Input", padding="10")
        csv_section.grid(row=0, column=0, sticky=(tk.W, tk.E), pady=(0, 10))
        csv_section.columnconfigure(0, weight=1)
        
        csv_input_frame = ttk.Frame(csv_section)
        csv_input_frame.grid(row=0, column=0, sticky=(tk.W, tk.E))
        csv_input_frame.columnconfigure(0, weight=1)
        
        self.csv_entry = ttk.Entry(csv_input_frame, textvariable=self.csv_file_path)
        self.csv_entry.grid(row=0, column=0, sticky=(tk.W, tk.E), padx=(0, 5))
        
        csv_browse_btn = ttk.Button(csv_input_frame, text="Browse CSV", 
                                    command=self.browse_csv)
        csv_browse_btn.grid(row=0, column=1)
        
        # CSV Preview
        csv_preview_label = ttk.Label(csv_section, text="Preview:", font=('Arial', 9, 'bold'))
        csv_preview_label.grid(row=1, column=0, sticky=tk.W, pady=(10, 2))
        
        csv_preview_frame = ttk.Frame(csv_section)
        csv_preview_frame.grid(row=2, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        csv_preview_frame.columnconfigure(0, weight=1)
        csv_preview_frame.rowconfigure(0, weight=1)
        
        self.csv_preview = tk.Text(csv_preview_frame, height=8, wrap=tk.NONE, 
                                   font=('Courier', 9), bg='#f5f5f5')
        self.csv_preview.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        
        csv_scroll_y = ttk.Scrollbar(csv_preview_frame, orient=tk.VERTICAL, 
                                     command=self.csv_preview.yview)
        csv_scroll_y.grid(row=0, column=1, sticky=(tk.N, tk.S))
        self.csv_preview.configure(yscrollcommand=csv_scroll_y.set)
        
        csv_scroll_x = ttk.Scrollbar(csv_preview_frame, orient=tk.HORIZONTAL, 
                                     command=self.csv_preview.xview)
        csv_scroll_x.grid(row=1, column=0, sticky=(tk.W, tk.E))
        self.csv_preview.configure(xscrollcommand=csv_scroll_x.set)
        
        csv_section.rowconfigure(2, weight=1)
        
        # TXT File Input Section
        txt_section = ttk.LabelFrame(left_frame, text="TXT File Input", padding="10")
        txt_section.grid(row=1, column=0, sticky=(tk.W, tk.E), pady=(0, 10))
        txt_section.columnconfigure(0, weight=1)
        
        txt_input_frame = ttk.Frame(txt_section)
        txt_input_frame.grid(row=0, column=0, sticky=(tk.W, tk.E))
        txt_input_frame.columnconfigure(0, weight=1)
        
        self.txt_entry = ttk.Entry(txt_input_frame, textvariable=self.txt_file_path)
        self.txt_entry.grid(row=0, column=0, sticky=(tk.W, tk.E), padx=(0, 5))
        
        txt_browse_btn = ttk.Button(txt_input_frame, text="Browse TXT", 
                                    command=self.browse_txt)
        txt_browse_btn.grid(row=0, column=1)
        
        # TXT Preview
        txt_preview_label = ttk.Label(txt_section, text="Preview:", font=('Arial', 9, 'bold'))
        txt_preview_label.grid(row=1, column=0, sticky=tk.W, pady=(10, 2))
        
        txt_preview_frame = ttk.Frame(txt_section)
        txt_preview_frame.grid(row=2, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        txt_preview_frame.columnconfigure(0, weight=1)
        txt_preview_frame.rowconfigure(0, weight=1)
        
        self.txt_preview = tk.Text(txt_preview_frame, height=8, wrap=tk.WORD, 
                                   font=('Courier', 9), bg='#f5f5f5')
        self.txt_preview.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        
        txt_scroll = ttk.Scrollbar(txt_preview_frame, orient=tk.VERTICAL, 
                                   command=self.txt_preview.yview)
        txt_scroll.grid(row=0, column=1, sticky=(tk.N, tk.S))
        self.txt_preview.configure(yscrollcommand=txt_scroll.set)
        
        txt_section.rowconfigure(2, weight=1)
        
        # Action Buttons
        button_frame = ttk.Frame(left_frame)
        button_frame.grid(row=2, column=0, pady=15)
        
        self.process_btn = ttk.Button(button_frame, text="Process Files", 
                                      command=self.process_files, width=18)
        self.process_btn.grid(row=0, column=0, padx=5)
        
        self.analyze_btn = ttk.Button(button_frame, text="Analyze Data", 
                                      command=self.analyze_data, width=18)
        self.analyze_btn.grid(row=0, column=1, padx=5)
        
        self.clear_btn = ttk.Button(button_frame, text="Clear All", 
                                    command=self.clear_all, width=18)
        self.clear_btn.grid(row=0, column=2, padx=5)
        
        # Table with 4 Options Section
        table_section = ttk.LabelFrame(left_frame, text="Output Options & Data", padding="10")
        table_section.grid(row=3, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        table_section.columnconfigure(0, weight=1)
        table_section.rowconfigure(1, weight=1)
        
        # Create notebook for tabs
        self.notebook = ttk.Notebook(table_section)
        self.notebook.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        
        # Create 4 option tabs
        self.option_frames = {}
        self.option_listboxes = {}
        # "fused_telomere","not_fused_telomere_(normal)","not_fused_no_telomere","fused_no_telomere"
        for i, option_name in enumerate(["Fused Telomere", "Not Fused Telomere (Normal)", "Not Fused No Telomere", "Fused No Telomere"], 1):
            tab_frame = ttk.Frame(self.notebook, padding="10")
            self.notebook.add(tab_frame, text=option_name)
            
            tab_frame.columnconfigure(0, weight=1)
            tab_frame.rowconfigure(0, weight=1)
            
            # Listbox for each option
            list_frame = ttk.Frame(tab_frame)
            list_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
            list_frame.columnconfigure(0, weight=1)
            list_frame.rowconfigure(0, weight=1)
            
            listbox = tk.Listbox(list_frame, font=('Courier', 9), height=10)
            listbox.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
            
            scrollbar = ttk.Scrollbar(list_frame, orient=tk.VERTICAL, command=listbox.yview)
            scrollbar.grid(row=0, column=1, sticky=(tk.N, tk.S))
            listbox.configure(yscrollcommand=scrollbar.set)
            
            self.option_listboxes[option_name] = listbox
            self.option_frames[option_name] = tab_frame
        
        # ===== RIGHT COLUMN =====
        right_frame = ttk.Frame(main_frame)
        right_frame.grid(row=1, column=1, rowspan=3, sticky=(tk.W, tk.E, tk.N, tk.S))
        right_frame.columnconfigure(0, weight=1)
        right_frame.rowconfigure(0, weight=1)
        
        # Image Output Section
        image_section = ttk.LabelFrame(right_frame, text="Output Alignment Visualized", padding="10")
        image_section.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        image_section.columnconfigure(0, weight=1)
        image_section.rowconfigure(0, weight=1)
        
        # Canvas for image display
        self.image_canvas = tk.Canvas(image_section, bg='white', 
                                      highlightthickness=1, highlightbackground='gray')
        self.image_canvas.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        
        # Image control buttons
        image_btn_frame = ttk.Frame(image_section)
        image_btn_frame.grid(row=1, column=0, pady=(10, 0))
        
        self.generate_img_btn = ttk.Button(image_btn_frame, text="Generate Image", 
                                          command=self.generate_output_image, width=18)
        self.generate_img_btn.grid(row=0, column=0, padx=5)
        
        self.save_img_btn = ttk.Button(image_btn_frame, text="Save Image", 
                                       command=self.save_image, width=18)
        self.save_img_btn.grid(row=0, column=1, padx=5)
        
        # Status bar
        self.status_var = tk.StringVar()
        self.status_var.set("Ready - Please select CSV and TXT files")
        status_bar = ttk.Label(self.root, textvariable=self.status_var, 
                              relief=tk.SUNKEN, anchor=tk.W)
        status_bar.grid(row=1, column=0, sticky=(tk.W, tk.E))
    
    def browse_csv(self):
        """Browse and select CSV file"""
        filename = filedialog.askopenfilename(
            title="Select a CSV file",
            filetypes=[("CSV files", "*.csv"), ("All files", "*.*")]
        )
        if filename:
            self.csv_file_path.set(filename)
            self.preview_csv(filename)
            self.update_status(f"CSV file selected: {os.path.basename(filename)}")
    
    def browse_txt(self):
        """Browse and select TXT file"""
        filename = filedialog.askopenfilename(
            title="Select a TXT file",
            filetypes=[("Text files", "*.txt"), ("All files", "*.*")]
        )
        if filename:
            self.txt_file_path.set(filename)
            self.preview_txt(filename)
            self.update_status(f"TXT file selected: {os.path.basename(filename)}")
    
    def preview_csv(self, filepath):
        """Preview CSV file content"""
        try:
            self.csv_preview.delete(1.0, tk.END)
            with open(filepath, 'r', encoding='utf-8') as f:
                reader = csv.reader(f)
                lines = []
                for i, row in enumerate(reader):
                    if i < 20:  # Show first 20 rows
                        lines.append(','.join(row))
                    else:
                        break
                
                preview_text = '\n'.join(lines)
                if i >= 20:
                    preview_text += '\n... (showing first 20 rows)'
                
                self.csv_preview.insert(1.0, preview_text)
        except Exception as e:
            self.csv_preview.insert(1.0, f"Error previewing CSV: {str(e)}")
    
    def preview_txt(self, filepath):
        """Preview TXT file content"""
        try:
            self.txt_preview.delete(1.0, tk.END)
            with open(filepath, 'r', encoding='utf-8') as f:
                content = f.read(5000)  # Read first 5000 characters
                self.txt_preview.insert(1.0, content)
                if len(content) >= 5000:
                    self.txt_preview.insert(tk.END, '\n... (showing first 5000 characters)')
        except Exception as e:
            self.txt_preview.insert(1.0, f"Error previewing TXT: {str(e)}")
    
    def process_files(self):
        """Process both files and populate option lists"""
        csv_file = self.csv_file_path.get()
        txt_file = self.txt_file_path.get()
        
        if not csv_file or not txt_file:
            messagebox.showwarning("Missing Files", 
                                  "Please select both CSV and TXT files!")
            return
        
        try:
            # Example processing - you'll replace this with your logic
            
            # Process CSV file
            csv_data = []
            with open(csv_file, 'r', encoding='utf-8') as f:
                reader = csv.reader(f)
                for row in reader:
                    csv_data.append(row)
            
            # Process TXT file
            with open(txt_file, 'r', encoding='utf-8') as f:
                txt_lines = f.readlines()
            
            # Populate Option 1: CSV headers or first row
            self.option_data["Option 1"] = csv_data[0] if csv_data else []
            self.update_option_list("Option 1")
            
            # Populate Option 2: CSV row count info
            self.option_data["Option 2"] = [
                f"Total CSV rows: {len(csv_data)}",
                f"Total columns: {len(csv_data[0]) if csv_data else 0}",
                f"First column: {csv_data[0][0] if csv_data and csv_data[0] else 'N/A'}"
            ]
            self.update_option_list("Option 2")
            
            # Populate Option 3: TXT file stats
            self.option_data["Option 3"] = [
                f"Total lines: {len(txt_lines)}",
                f"Total characters: {sum(len(line) for line in txt_lines)}",
                f"First line: {txt_lines[0].strip() if txt_lines else 'N/A'}"
            ]
            self.update_option_list("Option 3")
            
            # Populate Option 4: Combined info
            self.option_data["Option 4"] = [
                "=== Combined Summary ===",
                f"CSV file: {os.path.basename(csv_file)}",
                f"TXT file: {os.path.basename(txt_file)}",
                f"CSV rows: {len(csv_data)}",
                f"TXT lines: {len(txt_lines)}",
                "Status: Ready for analysis"
            ]
            self.update_option_list("Option 4")
            
            self.update_status("Files processed successfully")
            messagebox.showinfo("Success", "Files processed successfully!")
            
        except Exception as e:
            messagebox.showerror("Error", f"Error processing files: {str(e)}")
            self.update_status("Error processing files")
    
    def update_option_list(self, option_name):
        """Update the listbox for a specific option"""
        listbox = self.option_listboxes[option_name]
        listbox.delete(0, tk.END)
        for item in self.option_data[option_name]:
            listbox.insert(tk.END, item)
    
    def analyze_data(self):
        """
        Analyze data - this will tie into your external Python code
        PLACEHOLDER: Add your analysis code here
        """
        csv_file = self.csv_file_path.get()
        txt_file = self.txt_file_path.get()
        
        if not csv_file or not txt_file:
            messagebox.showwarning("Missing Files", 
                                  "Please select both files and process them first!")
            return
        
        try:
            # PLACEHOLDER: This is where you'll call your analysis code
            # Example:
            # from your_analysis_module import analyze_function
            # results = analyze_function(csv_file, txt_file)
            
            # For now, just show a message
            messagebox.showinfo("Analyze Data", 
                              "Analysis function will be integrated here.\n\n"
                              "This will call your external Python code with:\n"
                              f"CSV: {os.path.basename(csv_file)}\n"
                              f"TXT: {os.path.basename(txt_file)}")
            
            self.update_status("Analysis complete - ready to generate image")
            
        except Exception as e:
            messagebox.showerror("Error", f"Error analyzing data: {str(e)}")
    
    def generate_output_image(self):
        """Generate an image displaying all the output data"""
        try:
            # Create image
            img_width = 900
            img_height = 1000
            image = Image.new('RGB', (img_width, img_height), 'white')
            draw = ImageDraw.Draw(image)
            
            # Try to use a nicer font, fallback to default if not available
            try:
                title_font = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", 20)
                header_font = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", 14)
                text_font = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 11)
            except:
                title_font = ImageFont.load_default()
                header_font = ImageFont.load_default()
                text_font = ImageFont.load_default()
            
            # Title
            draw.text((20, 20), "File Processing Output Summary", fill='black', font=title_font)
            draw.line([(20, 50), (img_width-20, 50)], fill='gray', width=2)
            
            y_position = 70
            
            # File Information
            draw.text((20, y_position), "Input Files:", fill='blue', font=header_font)
            y_position += 25
            
            csv_name = os.path.basename(self.csv_file_path.get()) if self.csv_file_path.get() else "No CSV selected"
            txt_name = os.path.basename(self.txt_file_path.get()) if self.txt_file_path.get() else "No TXT selected"
            
            draw.text((40, y_position), f"CSV: {csv_name}", fill='black', font=text_font)
            y_position += 20
            draw.text((40, y_position), f"TXT: {txt_name}", fill='black', font=text_font)
            y_position += 40
            
            # Draw each option's data
            for option_name in ["Option 1", "Option 2", "Option 3", "Option 4"]:
                if y_position > img_height - 150:
                    break
                
                # Option header
                draw.text((20, y_position), f"{option_name}:", fill='blue', font=header_font)
                y_position += 25
                
                # Option data
                data_items = self.option_data.get(option_name, [])
                for item in data_items[:10]:  # Limit to 10 items per option
                    if y_position > img_height - 50:
                        break
                    item_str = str(item)[:100]  # Truncate long items
                    draw.text((40, y_position), f"• {item_str}", fill='black', font=text_font)
                    y_position += 18
                
                y_position += 20
            
            # Timestamp
            from datetime import datetime
            timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            draw.text((20, img_height - 30), f"Generated: {timestamp}", fill='gray', font=text_font)
            
            # Display image on canvas
            self.output_image = image
            self.display_image_on_canvas(image)
            
            self.update_status("Output image generated successfully")
            
        except Exception as e:
            messagebox.showerror("Error", f"Error generating image: {str(e)}")
            self.update_status("Error generating image")
    
    def display_image_on_canvas(self, image):
        """Display PIL image on canvas"""
        # Resize image to fit canvas if needed
        canvas_width = self.image_canvas.winfo_width()
        canvas_height = self.image_canvas.winfo_height()
        
        if canvas_width < 10:  # Canvas not yet rendered
            canvas_width = 800
            canvas_height = 900
        
        # Calculate scaling to fit canvas
        img_width, img_height = image.size
        scale = min(canvas_width / img_width, canvas_height / img_height, 1.0)
        
        new_width = int(img_width * scale)
        new_height = int(img_height * scale)
        
        resized_image = image.resize((new_width, new_height), Image.Resampling.LANCZOS)
        
        # Convert to PhotoImage
        self.photo = ImageTk.PhotoImage(resized_image)
        
        # Clear canvas and display image
        self.image_canvas.delete("all")
        self.image_canvas.create_image(canvas_width//2, canvas_height//2, 
                                       image=self.photo, anchor=tk.CENTER)
    
    def save_image(self):
        """Save the generated image"""
        if not self.output_image:
            messagebox.showwarning("No Image", "Please generate an image first!")
            return
        
        filename = filedialog.asksaveasfilename(
            title="Save Image As",
            defaultextension=".png",
            filetypes=[("PNG files", "*.png"), ("JPEG files", "*.jpg"), ("All files", "*.*")]
        )
        
        if filename:
            try:
                self.output_image.save(filename)
                messagebox.showinfo("Success", f"Image saved to:\n{filename}")
                self.update_status(f"Image saved: {os.path.basename(filename)}")
            except Exception as e:
                messagebox.showerror("Error", f"Error saving image: {str(e)}")
    
    def clear_all(self):
        """Clear all inputs and outputs"""
        self.csv_file_path.set("")
        self.txt_file_path.set("")
        self.csv_preview.delete(1.0, tk.END)
        self.txt_preview.delete(1.0, tk.END)
        
        for option_name in self.option_data:
            self.option_data[option_name] = []
            self.update_option_list(option_name)
        
        self.image_canvas.delete("all")
        self.output_image = None
        
        self.update_status("All data cleared")
    
    def update_status(self, message):
        """Update the status bar"""
        self.status_var.set(message)


def main():
    root = tk.Tk()
    app = EnhancedFileApp(root)
    root.mainloop()


if __name__ == "__main__":
    main()
