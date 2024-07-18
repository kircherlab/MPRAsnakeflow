from IPython.display import Image, display

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
        print(png_file)
        display(Image(png_file))
        # print(png_file)
