#!/bin/bash
filename=${0%.sh}
job_dir=`dirname "$(realpath $0)"`
docker_image="europe-southwest1-docker.pkg.dev/repositories-368911/proteinprepare-v1/proteinprepare:latest"
license_server_port=27000
exposes_ports=

# Function to check if a port is in use
is_port_in_use() {
    local port=$1
    
    # First check if Docker has the port allocated
    if command -v docker &> /dev/null; then
        docker ps --format "table {{.Ports}}" 2>/dev/null | grep -q "0.0.0.0:$port->" && return 0
        docker ps --format "table {{.Ports}}" 2>/dev/null | grep -q ":$port->" && return 0
    fi
    
    # Then check system ports
    if command -v lsof &> /dev/null; then
        lsof -i:$port -sTCP:LISTEN -t >/dev/null 2>&1
        return $?
    elif command -v nc &> /dev/null; then
        nc -z localhost $port >/dev/null 2>&1
        return $?
    elif command -v ss &> /dev/null; then
        ss -tuln | grep -q ":$port "
        return $?
    elif command -v netstat &> /dev/null; then
        netstat -tuln | grep -q ":$port "
        return $?
    else
        echo "Warning: No port checking utility found (lsof, nc, ss, netstat). Assuming port is available."
        return 1
    fi
}

# Find available ports and build port mapping flags
port_flags=""
port_env_flags=""
max_attempts=100

if [ -n "$exposes_ports" ]; then
    # Split the comma-separated ports
    IFS=',' read -ra PORTS <<< "$exposes_ports"
    
    # Track used host ports to avoid collisions when multiple container ports map to incremented host ports
    declare -A used_host_ports
    
    for container_port in "${PORTS[@]}"; do
        # Trim whitespace
        container_port=$(echo "$container_port" | xargs)
        
        # Skip empty values
        if [ -z "$container_port" ]; then
            continue
        fi
        
        # Find an available host port starting from the container port number
        host_port=$container_port
        attempts=0
        
        while is_port_in_use $host_port || [ -n "${used_host_ports[$host_port]}" ]; do
            if [ $attempts -ge $max_attempts ]; then
                echo "Error: Could not find an available port after $max_attempts attempts for container port $container_port"
                exit 1
            fi
            echo "Port $host_port is already in use, trying next port..."
            host_port=$((host_port + 1))
            attempts=$((attempts + 1))
        done
        
        # Mark this host port as used
        used_host_ports[$host_port]=1
        
        echo "Mapping container port $container_port to host port $host_port"
        
        # Add port mapping flag
        port_flags="$port_flags -p $host_port:$container_port"
        
        # Add environment variable for the port mapping
        # Format: PORT_<container_port>=<host_port>
        port_env_flags="$port_env_flags -e PORT_${container_port}=$host_port"
    done
fi

# flags needed for the license checker
docker_flags="-v /etc/passwd:/etc/passwd:ro \
    -v /etc/group:/etc/group:ro \
    --add-host=host.docker.internal:host-gateway \
    -e ACELLERA_LICENCE_SERVER=${ACELLERA_LICENCE_SERVER:-$license_server_port@host.docker.internal} \
    --shm-size=8g" 

# GPU support: handle CUDA_VISIBLE_DEVICES if set
gpu_flags=""
if [ -n "$CUDA_VISIBLE_DEVICES" ]; then
    # CUDA_VISIBLE_DEVICES is set, expose those specific GPUs (they will be renumbered from 0 inside the container)
    gpu_flags="--gpus \"device=$CUDA_VISIBLE_DEVICES\""
else
    # CUDA_VISIBLE_DEVICES not set, expose all GPUs
    gpu_flags="--gpus all"
fi

# If the CI envvar is set, pass it to the container
# Note: --init flag ensures proper signal handling and process reaping
# When this script is terminated (e.g., via Ctrl+C), the signal will be properly
# propagated to the container. The --rm flag ensures the container is automatically removed.
exec docker run --init --user $(id -u):$(id -g) --rm \
    $docker_flags $gpu_flags $port_flags $port_env_flags \
    -v $job_dir:/workspace \
    $docker_image \
    --input-json $filename/inputs.json /workspace

exit_code=$?

if [ $exit_code -ne 0 ]; then
    echo "Error in running application. Check logs..."
    exit $exit_code
fi