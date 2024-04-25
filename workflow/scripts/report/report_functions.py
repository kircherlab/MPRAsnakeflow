from IPython.display import Image, display

def get_unique_files(png_files):
    file_dict = {}

    # Loop through the files
    for file in png_files:
        # If the file name is not in the dictionary, add it
        if file.name not in file_dict:
              file_dict[file.name] = file

    # Get the unique files
    unique_files = list(file_dict.values())
    return unique_files


def display_unique_images(pattern, image_path, exclude=None):
    # Get the image files that match the pattern
    png_files = list(image_path.glob(pattern))
    
    # Exclude files if necessary
    if exclude:
        png_files = [file for file in png_files if exclude not in file.name]

    # Get the unique files
    unique_files = get_unique_files(png_files)

    # Display the images and print their paths
    for png_file in unique_files:
        display(Image(png_file))
        # print(png_file)