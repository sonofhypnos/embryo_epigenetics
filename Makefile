# Define the Conda environment name
CONDA_ENV := epi_env

# Specify the Python script
SCRIPT := src/cosine_similarity.py

# Path to the results directory
OUTPUT_DIR := results

# Default rule that runs when you simply type "make" in the terminal
all: run_script

# Rule to create the output directory if it does not exist
$(OUTPUT_DIR):
	mkdir -p $@

# Rule to run the Python script within the Conda environment
run_script: | $(OUTPUT_DIR)
	@echo "Activating Conda environment $(CONDA_ENV) and running script..."
	@SHELL=/bin/zsh conda run -n $(CONDA_ENV) python $(SCRIPT)

# Clean rule to potentially remove the results directory
clean:
	rm -rf $(OUTPUT_DIR)

.PHONY: all run_script clean
