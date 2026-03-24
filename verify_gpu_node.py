import torch
import os

def verify_gpu():
    print("=== GPU Assignment Verification ===")
    print(f"Process ID: {os.getpid()}")
    print(f"CUDA_VISIBLE_DEVICES: {os.environ.get('CUDA_VISIBLE_DEVICES', 'Not Set')}")
    
    if not torch.cuda.is_available():
        print("❌ CUDA is not available to this process.")
        return

    device_count = torch.cuda.device_count()
    current_device = torch.cuda.current_device()
    device_name = torch.cuda.get_device_name(current_device)
    
    print(f"📝 CUDA Device Count: {device_count}")
    print(f"✅ Currently Using: Device {current_device} ({device_name})")
    print(f"💾 VRAM Allocated: {torch.cuda.memory_allocated(current_device) / 1024**2:.2f} MB")
    print("====================================")

if __name__ == "__main__":
    verify_gpu()
