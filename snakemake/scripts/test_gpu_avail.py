### CHECK GPU AVAILABILITY ###
import torch
import GPUtil
print('GPU available:', GPUtil.getAvailable())
print('CUDA available:', torch.cuda.is_available())
if len(GPUtil.getAvailable()) > 0: print('CUDA device:', torch.cuda.get_device_name(0))