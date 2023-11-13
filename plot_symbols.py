import matplotlib.pyplot as plt
import numpy as np
import cv2

# Example usage:
image_freezer = cv2.imread("freezer.png")  # Load your image here
image_heater = cv2.imread("heater.png")  # Load your image here
image = cv2.imread("normal_transfer.png")  # Load your image here


text_labels = {
    (0.55 * image.shape[1], 0.05 * image.shape[0]): "1",
    (0.05 * image.shape[1], 0.55 * image.shape[0]): "2",
    (0.95 * image.shape[1], 0.55 * image.shape[0]): "3",
    (0.55 * image.shape[1], 0.95 * image.shape[0]): "4",
    (0.5 * image.shape[1], 0.5 * image.shape[0]): "1"
}

text_labels_heater = {
    (0.55 * image.shape[1], 0.05 * image.shape[0]): "1",
    (0.25 * image.shape[1], 0.8 * image.shape[0]): "2",
    (0.75 * image.shape[1], 0.2 * image.shape[0]): "3",
    (0.55 * image.shape[1], 0.95 * image.shape[0]): "4",
    (0.5 * image.shape[1], 0.5 * image.shape[0]): "1"
}

text_labels_freezer = {
    (0.05 * image.shape[1], 0.57 * image.shape[0]): "1",
    (0.25 * image.shape[1], 0.8 * image.shape[0]): "2",
    (0.8 * image.shape[1], 0.2 * image.shape[0]): "3",
    (0.95 * image.shape[1], 0.55 * image.shape[0]): "4",
    (0.5 * image.shape[1], 0.5 * image.shape[0]): "1"
}

def plot_image_with_numbers(image, text_list):
    """
    Plot an image and display text next to each section as shown in the diagram.

    Args:
    image (numpy.ndarray): The input image to be displayed.
    text_list (dict): A dictionary where keys are position tuples and values are the text labels.

    Returns:
    None
    """
    fig, ax = plt.subplots()

    # Display the image
    ax.imshow(image, cmap='gray')

    for pos, text in text_list.items():
        ax.text(pos[0], pos[1], text, fontsize=12, color='red',
                ha='center', va='center')

    # Remove axis labels and ticks
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    
    return image
    #plt.show()

import numpy as np
import matplotlib.pyplot as plt

def combine_images(images, orientation='horizontal'):
    if not images:
        return None

    if orientation not in ['horizontal', 'vertical']:
        raise ValueError("Invalid orientation. Use 'horizontal' or 'vertical'.")

    if orientation == 'horizontal':
        total_height = max(image.shape[0] for image in images)
        total_width = sum(image.shape[1] for image in images)
        combined_image = np.zeros((total_height, total_width, 3), dtype=np.uint8)

        x_offset = 0
        for image in images:
            if image.shape[2] == 4:  # Check for an RGBA image
                # If the image has 4 color channels (RGBA), remove the alpha channel
                image = image[:, :, :3]
            combined_image[:image.shape[0], x_offset:x_offset + image.shape[1], :] = image
            x_offset += image.shape[1]

    elif orientation == 'vertical':
        total_width = max(image.shape[1] for image in images)
        total_height = sum(image.shape[0] for image in images)
        combined_image = np.zeros((total_height, total_width, 3), dtype=np.uint8)

        y_offset = 0
        for image in images:
            if image.shape[2] == 4:
                image = image[:, :, :3]
            combined_image[y_offset:y_offset + image.shape[0], :image.shape[1], :] = image
            y_offset += image.shape[0]

    # Swap the red and blue channels
    combined_image[:, :, [0, 2]] = combined_image[:, :, [2, 0]]

    plt.imshow(combined_image)
    plt.axis('off')
    plt.show()

combine_images([image_heater, image], orientation='vertical')

#plot_image_with_numbers(image_freezer, text_labels_freezer)
