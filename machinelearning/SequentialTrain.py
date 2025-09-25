import torch
import torch.nn as nn
import torch.optim as optim
import torch.distributed as dist
import torch.multiprocessing as mp
import os

# --- Configuration ---
WORLD_SIZE = 4
MASTER_PORT = '29500' 
DATASET_SIZE_PER_RANK = 1000
INPUT_FEATURES = 10
LEARNING_RATE = 0.01
NUM_EPOCHS = 5
BATCH_SIZE = 16
# ---------------------

# --- 1. Define the Simple Model ---
class SimpleNN(nn.Module):
    def __init__(self):
        super(SimpleNN, self).__init__()
        self.layer1 = nn.Linear(INPUT_FEATURES, 20)
        self.relu = nn.ReLU()
        self.layer2 = nn.Linear(20, 1)

    def forward(self, x):
        return self.layer2(self.relu(self.layer1(x)))

# --- PHASE 1 Worker: Dataset Creation ---
def dataset_creator_worker(rank, world_size, data_tensor, target_tensor):
    """
    Worker function to create a unique chunk of the shared dataset.
    """
    print(f"[{rank}/{world_size}] PHASE 1: Creating data chunk...")
    
    # Calculate the starting index for this rank's chunk
    start_idx = rank * DATASET_SIZE_PER_RANK
    end_idx = start_idx + DATASET_SIZE_PER_RANK

    # Generate mock data
    # Input data (add rank offset to ensure uniqueness)
    data_chunk = torch.randn(DATASET_SIZE_PER_RANK, INPUT_FEATURES) + rank * 0.1
    # Target data (sum of inputs)
    target_chunk = data_chunk.sum(dim=1, keepdim=True)

    # Place chunks into the shared tensors
    data_tensor[start_idx:end_idx] = data_chunk
    target_tensor[start_idx:end_idx] = target_chunk
    
    print(f"[{rank}/{world_size}] PHASE 1: Data chunk saved to shared memory.")

# --- PHASE 2 Worker: DDP Training ---
def ddp_training_worker(rank, world_size, data_tensor, target_tensor):
    """
    Worker function to initialize DDP and train the model using the shared dataset.
    """
    
    # 1. Setup Environment and Initialize Process Group
    os.environ['MASTER_ADDR'] = 'localhost'
    os.environ['MASTER_PORT'] = MASTER_PORT
    
    dist.init_process_group(
        backend="gloo",
        rank=rank,
        world_size=world_size
    )
    print(f"[{rank}/{world_size}] PHASE 2: DDP group initialized.")

    # 2. Setup Distributed Data Sampler
    # We use a standard TensorDataset and DistributedSampler for DDP compatibility
    full_dataset = torch.utils.data.TensorDataset(data_tensor, target_tensor)
    sampler = torch.utils.data.distributed.DistributedSampler(
        full_dataset, 
        num_replicas=world_size, 
        rank=rank, 
        shuffle=True
    )
    dataloader = torch.utils.data.DataLoader(
        full_dataset, 
        batch_size=BATCH_SIZE, 
        sampler=sampler
    )
    
    # 3. Setup Model, Optimizer, and Loss
    model = SimpleNN()
    ddp_model = nn.parallel.DistributedDataParallel(model)
    criterion = nn.MSELoss()
    optimizer = optim.SGD(ddp_model.parameters(), lr=LEARNING_RATE)

    # 4. Training Loop
    for epoch in range(NUM_EPOCHS):
        # Must set sampler epoch so data is shuffled differently each epoch
        sampler.set_epoch(epoch)
        
        for batch_idx, (data, target) in enumerate(dataloader):
            optimizer.zero_grad()
            output = ddp_model(data)
            loss = criterion(output, target)
            loss.backward()
            optimizer.step()

        if rank == 0:
             print(f"[Rank {rank}/{world_size}] Epoch {epoch+1}/{NUM_EPOCHS} Loss: {loss.item():.4f}")

    # 5. Clean up
    dist.destroy_process_group()


if __name__ == '__main__':
    # Total size of the dataset
    TOTAL_DATASET_SIZE = WORLD_SIZE * DATASET_SIZE_PER_RANK
    
    # Set the start method for robustness
    try:
        mp.set_start_method('spawn', force=True)
    except RuntimeError:
        pass 

    # --- Setup Shared Memory Tensors (Global Dataset) ---
    data_shared = torch.zeros(TOTAL_DATASET_SIZE, INPUT_FEATURES).share_memory_()
    target_shared = torch.zeros(TOTAL_DATASET_SIZE, 1).share_memory_()

    # =================================================================
    # PHASE 1: CREATE SHARED DATASET (Multiprocessing)
    # =================================================================
    print("\n--- Starting Phase 1: Creating Shared Dataset ---")
    mp.spawn(
        dataset_creator_worker,
        args=(WORLD_SIZE, data_shared, target_shared),
        nprocs=WORLD_SIZE,
        join=True  # Blocks until all data workers have terminated
    )
    print("\n--- Phase 1 Complete: Shared Dataset is Ready ---")
    
    # =================================================================
    # PHASE 2: DDP TRAINING (Distributed)
    # =================================================================
    print("\n--- Starting Phase 2: DDP Training ---")
    mp.spawn(
        ddp_training_worker,
        args=(WORLD_SIZE, data_shared, target_shared),
        nprocs=WORLD_SIZE,
        join=True  # Blocks until all DDP workers have terminated
    )

    print("\n--- All Execution Complete ---")
    # The data_shared and target_shared tensors now hold the same data 
    # but the model is trained and its weights are gone when processes terminate.
