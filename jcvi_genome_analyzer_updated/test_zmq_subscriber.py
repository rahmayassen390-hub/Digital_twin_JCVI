#!/usr/bin/env python3
"""
Test ZeroMQ Subscriber for JCVI Genome Analyzer
Run this script to test if ZMQ broadcasting is working

Usage: python test_zmq_subscriber.py
"""

try:
    import zmq
    import json
except ImportError:
    print("Error: pyzmq not installed")
    print("Install with: pip install pyzmq")
    exit(1)

def main():
    context = zmq.Context()
    socket = context.socket(zmq.SUB)
    socket.connect("tcp://localhost:5555")
    socket.setsockopt_string(zmq.SUBSCRIBE, "")  # Subscribe to all messages
    
    print("🟢 ZMQ Subscriber connected to localhost:5555")
    print("Waiting for translation data...\n")
    print("=" * 70)
    
    data_count = 0
    while True:
        try:
            message = socket.recv_string()
            msg = json.loads(message)
            
            msg_type = msg.get('type', 'data')
            
            if msg_type == 'status':
                status = msg.get('status', 'unknown')
                timestamp = msg.get('timestamp', 0)
                print(f"\n{'='*70}")
                print(f"📢 STATUS: {status.upper()}")
                print(f"   Timestamp: {timestamp:.1f}s")
                print(f"{'='*70}\n")
                
                if status == 'completed':
                    print("✅ Animation completed - data is now stable")
                elif status == 'started':
                    print("🚀 Animation started - data streaming...")
                    data_count = 0
                elif status == 'reset':
                    print("🔄 Reset - clearing data")
                    data_count = 0
            
            else:
                data_count += 1
                print(f"[{data_count}] Gene: {msg['gene_id']} | Status: {msg['status']}")
                
                # Transcription info
                mrna_cur = msg.get('mrna_current', 0)
                mrna_tot = msg.get('mrna_copies', 0)
                tx_rate = msg.get('transcription_rate', 0)
                print(f"    📝 Transcription: mRNA {mrna_cur:.0f}/{mrna_tot} at {tx_rate:.2f}/sec")
                
                # Translation info
                prot_cur = msg.get('protein_current', 0)
                prot_tot = msg.get('protein_total', 0)
                tl_rate = msg.get('translation_rate', 0)
                print(f"    🧪 Translation:  Proteins {prot_cur:.0f}/{prot_tot} at {tl_rate:.1f}/sec")
                
                # Gene info
                gene_len = msg.get('gene_length', 0)
                prot_len = msg.get('protein_length', 0)
                promoter = msg.get('promoter_type', 'None')
                print(f"    📋 Gene: {gene_len}bp → {prot_len}aa | Promoter: {promoter}")
                print("-" * 50)
            
        except KeyboardInterrupt:
            print("\n\n🔴 Subscriber stopped")
            break
    
    socket.close()
    context.term()

if __name__ == "__main__":
    main()
