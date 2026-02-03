from moviepy.editor import *

image_folder = ''  
output_video = 'video.mp4'  
fps = 15  # Frames per second
image_files = [f"{image_folder}{str(i).zfill(3)}.jpg" for i in range(300)]
clip = ImageSequenceClip(image_files, fps=fps)

# write the video
clip.write_videofile(output_video, codec='libx264')
