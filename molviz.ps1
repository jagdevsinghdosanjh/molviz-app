# setup.ps1 — Initialize molviz-app environment

Write-Host "Setting up molviz-app..."

# Create virtual environment
python -m venv .venv
Write-Host "Virtual environment created."

# Activate virtual environment
& .\.venv\Scripts\Activate.ps1
Write-Host "Virtual environment activated."

# Install dependencies
pip install -r requirements.txt
Write-Host "Dependencies installed."

# Create folders
New-Item -ItemType Directory -Path "app"
New-Item -ItemType Directory -Path "assets"
New-Item -ItemType Directory -Path "assets/icons"
New-Item -ItemType Directory -Path "tests"
Write-Host "Project folders scaffolded."

# Create placeholder files
New-Item -ItemType File -Path "app/__init__.py"
New-Item -ItemType File -Path "app/main.py"
New-Item -ItemType File -Path "README.md"
New-Item -ItemType File -Path ".gitignore"
Write-Host "Starter files created."

Write-Host "Setup complete. You’re ready to build!"
