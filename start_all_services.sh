#!/bin/bash

# OQTOPUS Master Startup Script
# Starts Backend (Docker), APIs, Frontend, and Mock Worker

echo "🐙 Starting OQTOPUS Quantum Environment..."

# 1. Start Backend Services (Docker)
echo "\\n🐳 Starting Backend Services (MySQL, MinIO)..."
cd /home/rahema/Desktop/oqtopus-cloud/backend

# Fix for "ContainerConfig" KeyErrors in older docker-compose
echo "🧹 Cleaning up old containers to prevent configuration errors..."
docker-compose down --remove-orphans > /dev/null 2>&1
# Force remove if down failed to clean up partially created containers
docker rm -f backend_db_1 backend_minio_1 backend_minio-client_1 2>/dev/null

echo "🚀 Bringing up Docker services..."
docker-compose up -d
if [ $? -eq 0 ]; then
    echo "✅ Backend Services Started."
else
    echo "❌ Failed to start Backend Services. Retrying with --force-recreate..."
    docker-compose up -d --force-recreate
    if [ $? -eq 0 ]; then
        echo "✅ Backend Services Started (after retry)."
    else
        echo "❌ Critical Error: Could not start Backend Services."
        exit 1
    fi
fi

# Environment variables for local development
export ENV=local
export SECRET_NAME=local
export DB_HOST=localhost
export DB_NAME=main
export DB_CONNECTOR=mysql+pymysql
export MYSQL_USER=admin
export MYSQL_PASSWORD=password

# 2. Start User API (Port 8081)
echo "\\n🔌 Check/Start User API (Port 8081)..."
if lsof -i :8081 > /dev/null; then
    echo "✅ User API is already running."
else
    echo "🚀 Starting User API..."
    # Running in background with env vars
    nohup bash -c 'cd /home/rahema/Desktop/oqtopus-cloud/backend && ENV=local SECRET_NAME=local DB_HOST=localhost DB_NAME=main DB_CONNECTOR=mysql+pymysql MYSQL_USER=admin MYSQL_PASSWORD=password poetry run uvicorn oqtopus_cloud.user.lambda_function:app --host 0.0.0.0 --port 8081' > /tmp/oqtopus_user_api.log 2>&1 &
    echo "✅ User API started in background."
fi

# 3. Start Provider API (Port 8888)
echo "\\n🔌 Check/Start Provider API (Port 8888)..."
if lsof -i :8888 > /dev/null; then
    echo "✅ Provider API is already running."
else
    echo "🚀 Starting Provider API..."
    nohup bash -c 'cd /home/rahema/Desktop/oqtopus-cloud/backend && ENV=local SECRET_NAME=local DB_HOST=localhost DB_NAME=main DB_CONNECTOR=mysql+pymysql MYSQL_USER=admin MYSQL_PASSWORD=password poetry run uvicorn oqtopus_cloud.provider.lambda_function:app --host 0.0.0.0 --port 8888' > /tmp/oqtopus_provider_api.log 2>&1 &
    echo "✅ Provider API started in background."
fi

# 4. Start Frontend (Port 5173)
echo "\\n💻 Check/Start Frontend (Port 5173)..."
if lsof -i :5173 > /dev/null; then
    echo "✅ Frontend is already running."
else
    echo "🚀 Starting Frontend..."
    cd /home/rahema/Desktop/oqtopus-frontend
    nohup npm run dev > /tmp/oqtopus_frontend.log 2>&1 &
    echo "✅ Frontend started in background."
fi

# 5. Start Mock Simulator Worker
echo "\\n🤖 Check/Start Mock Simulator Worker..."
# Simple check if python process with mock_simulator_worker.py is running
if pgrep -f "mock_simulator_worker.py" > /dev/null; then
    echo "✅ Mock Simulator Worker is already running."
else
    echo "🚀 Starting Mock Simulator Worker..."
    cd "/home/rahema/Desktop/bed_files /jcvi_genome_analyzer_updated (6)"
    nohup python3 mock_simulator_worker.py > /tmp/oqtopus_mock_worker.log 2>&1 &
    echo "✅ Mock Simulator Worker started in background."
fi

echo "\\n🐙 All OQTOPUS services checked/started!"
echo "User API Logs: /tmp/oqtopus_user_api.log"
echo "Provider API Logs: /tmp/oqtopus_provider_api.log"
echo "Frontend Logs: /tmp/oqtopus_frontend.log"
echo "Mock Worker Logs: /tmp/oqtopus_mock_worker.log"
