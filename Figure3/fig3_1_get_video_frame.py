# Install Dependencies
import os
import cv2
import matplotlib.pyplot as plt

# Class Object for Video Frame Extraction
class VideoDecoder():
    def __init__(self, file_path, input_time):
        self.file_path = file_path
        self.input_time = input_time
        
    def extract_video_instance(self, save=False, save_path=None):
        videocap = cv2.VideoCapture(self.file_path)
        FPS = int(videocap.get(cv2.CAP_PROP_FPS))
        TOTAL_FRAME = int(videocap.get(cv2.CAP_PROP_FRAME_COUNT))
        CAP_PROP_POS_FRAMES = 1 # use 0-based index for decoding the frames
        frame_idx = self.input_time * FPS
        videocap.set(CAP_PROP_POS_FRAMES, frame_idx) # captures `frame_idx + 1`-th frame
        ret, frame = videocap.read()
        if not ret:
            raise ValueError('Return value is false.')
        print("Extracting Video Frame ... ")
        print("FPS: {} | Total frames: {} | Current Time: {} s | Current Frame: {}".format(FPS, TOTAL_FRAME, self.input_time, frame_idx))
        if save:
            if save_path is None:
                save_path = os.path.basename(self.file_path)
            cv2.imwrite(os.path.join(save_path, 'example_frame.png'), frame)
        else:
            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 6))
            plt.imshow(frame)
            ax.axis('off')
            plt.show()
        print("Extraction Complete.")
        return None

if __name__ == "__main__":
    # [1] Set Input Parameters
    input_time = 85.9326 # time (s) in video to be captured
    # NOTE: The times in signal and video were synchronized based on their respective recording 'start' times.

    # [2] Set System Parameters
    file_path = os.path.join('/Users/scho/Neuroscience_KIST/Cho2022_BurstDetection'
                            '/data/experimental_data/ESCAPE/sample_videos', 'escape_811.mp4')
    save_path = os.getcwd()
    if not os.path.exists(file_path):
        raise ValueError('Your video file does not exist.')

    # [3] Extract Video Frame
    VD = VideoDecoder(file_path, input_time)
    VD.extract_video_instance(save=True, save_path=save_path)