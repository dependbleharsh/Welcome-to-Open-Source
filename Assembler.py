import re
import sys

import pandas as pd

instructions_dict = {
    # J-type
    "jal": ("1101111", None, None),
    # U-type
    "lui": ("0110111", None, None),
    "auipc": ("0010111", None, None),
    # B-type
    "beq": ("1100011", "000", None),
    "bne": ("1100011", "001", None),
    "blt": ("1100011", "100", None),
    "bge": ("1100011", "101", None),
    "bltu": ("1100011", "110", None),
    "bgeu": ("1100011", "111", None),
    # S-type
    "sw": ("0100011", "010", None),
    # I-type
    "lw": ("0000011", "010", None),
    "addi": ("0010011", "000", None),
    "sltiu": ("0010011", "011", None),
    "jalr": ("1100111", "000", None),
    # R-type
    "add": ("0110011", "000", "0000000"),
    "sub": ("0110011", "000", "0100000"),
    "sll": ("0110011", "001", "0000000"),
    "slt": ("0110011", "010", "0000000"),
    "sltu": ("0110011", "011", "0000000"),
    "xor": ("0110011", "100", "0000000"),
    "srl": ("0110011", "101", "0000000"),
    "or": ("0110011", "110", "0000000"),
    "and": ("0110011", "111", "0000000")
}

registers = {
    "zero": "00000", "ra": "00001", "sp": "00010", "gp": "00011", "tp": "00100",
    "t0": "00101", "t1": "00110", "t2": "00111", "s0": "01000", "s1": "01001",
    "a0": "01010", "a1": "01011", "a2": "01100", "a3": "01101", "a4": "01110",
    "a5": "01111", "a6": "10000", "a7": "10001", "s2": "10010", "s3": "10011",
    "s4": "10100", "s5": "10101", "s6": "10110", "s7": "10111", "s8": "11000",
    "s9": "11001", "s10": "11010", "s11": "11011", "t3": "11100", "t4": "11101",
    "t5": "11110", "t6": "11111"
}


class Assembler:
    def __init__(self, filename):
        self.instructions = []
        with open(filename, "r") as f:
            self.lines = f.readlines()
            self.instructions = [line.strip() for line in self.lines if line]
        print(self.instructions)
        
        self.labels=[]
        for i in range(len(self.instructions)):
            if(re.search('\:',self.instructions[i])):
                split=re.split(":",self.instructions[i])
                #check for trailing white spaces in label
                if(split[0].strip()!=split[0]):
                    print("Invalid label at address ", 4*i)
                    exit()
                self.labels.append([split[0], 4*i])
                self.instructions[i] = split[1].strip()
        print(self.labels)
                
        self.asm1 = []
        for i in range(len(self.instructions)):
            self.asm1.append(re.split(" ", self.instructions[i], 1))
            
        asmdf = pd.DataFrame(self.asm1)
        print(asmdf)
        
        #search for "beq zero,zero,0" in the instructions
        a=self.asm1.index(["beq","zero,zero,0"])
        print("Index of virtual halt instruction: ",a)      #changed that virtual halt need not be the last instruction
        if(a==-1):
            print("Virtual halt instruction not found")
            #return
            exit()
        
    def assemble_instruction(self, instruction, current_address):
        # if(instruction == None): 
        #     return
        if(instruction[0] not in instructions_dict):
            print("Invalid instruction at address ", current_address)
            #return
            exit()
        opcode, funct3, funct7 = instructions_dict[instruction[0]]
        arguments = [arg.strip() for arg in re.split(",|\(|\)|:", instruction[1]) if arg]
        print(arguments)
        if(len(instruction)!=2):
            print("Invalid instruction at address ", current_address)
            #return
            exit()

        if instruction[0] in ["lui", "auipc"]:
            if(len(arguments)!=2):
                print("Invalid instruction at address ", current_address)
                #return
                exit()
            
            if(arguments[0] not in registers):
                print("Invalid register at address ", current_address)
                exit()
                #return
            rd = registers[arguments[0]]
            
            if(int(arguments[1]) >= 2**31 or int(arguments[1]) < -2**31):
                print("Invalid immediate value")
                exit()
                #return
            
            if(int(arguments[1]) >= 0):
                imm = '{:032b}'.format(int(arguments[1]))
            else:
                imm = '{:032b}'.format(2**32 + int(arguments[1]))
            return imm[0:20] + str(rd) + str(opcode)
        
        if instruction[0] in ["jal"]:
            if(len(arguments)!=2):
                print("Invalid instruction at address ", current_address)
                #return
                exit()
            
            if(arguments[0] not in registers):
                print("Invalid register at address ", current_address)
                exit()
            rd = registers[arguments[0]]
            try:        #check for immediate value
                imm = int(arguments[1])
                if(imm >= 2**20 or imm < -2**20):
                    print("Invalid immediate value")
                    return
                    exit()
                if(imm >= 0):
                    imm = '{:021b}'.format(imm)
                else:
                    imm = '{:021b}'.format(2**21 + imm)
            except:     #else check for label
                for i in self.labels:
                    if(i[0] == arguments[1]):
                        if(i[1]>=self.current_address):
                            imm = '{:021b}'.format(int((i[1] - current_address)))
                        else:
                            imm = '{:021b}'.format(2**21 + int((i[1] - current_address)))
                        break
                else:    
                    print("Invalid label at address ", current_address)
            return imm[0] + imm[-11:-1] + imm[10] + imm[1:9] + str(rd) + str(opcode)

        if instruction[0] in ["jalr"]:
            if(len(arguments)!=3):
                print("Invalid instruction at address ", current_address)
                exit()
            
            if(arguments[0] not in registers or arguments[1] not in registers):
                print("Invalid register at address ", current_address)
                exit()
            
            rd = registers[arguments[0]]
            rs1 = registers[arguments[1]]
            
            try:        #check for immediate value
                imm = int(arguments[2])
                if(imm >= 2**11 or imm < -2**11):
                    print("Invalid immediate value")
                    return
                    exit()
                if(imm >= 0):
                    imm = '{:012b}'.format(imm)
                else:
                    imm = '{:012b}'.format(2**12 + imm)
            except:     #else check for label
                for i in self.labels:
                    if(i[0] == arguments[2]):
                        if(i[1]>=self.current_address):
                            imm = '{:012b}'.format(int((i[1] - current_address)))
                        else:
                            imm = '{:012b}'.format(2**12 + int((i[1] - current_address)))
                        break
                else:
                    print("Invalid label at address ", current_address)
                    exit()
                
            return imm + str(rs1) + str(funct3) + str(rd) + str(opcode)

        if instruction[0] in ["beq", "bne", "blt", "bge", "bltu", "bgeu"]:
            if(len(arguments)!=3):
                print("Invalid instruction at address ", current_address)
                exit()
            
            if(arguments[0] not in registers or arguments[1] not in registers):
                print("Invalid register at address ", current_address)
                exit()
            
            rs1 = registers[arguments[0]]
            rs2 = registers[arguments[1]]
            
            try:            #check for immediate value
                imm = int(arguments[2])
                if(imm >= 2**12 or imm < -2**12):
                    print("Invalid immediate value")
                    return
                    exit()
                    #return
                if(imm >= 0):
                    imm = '{:013b}'.format(imm)
                else:
                    imm = '{:013b}'.format(2**13 + imm)
            except:        #else check for label
                for i in self.labels:
                    if(i[0] == arguments[2]):
                        if(i[1]>=self.current_address):
                            imm = '{:013b}'.format(int((i[1] - current_address)))
                        else:
                            imm = '{:013b}'.format(2**13 + int((i[1] - current_address)))
                        break
                else:
                    print("Invalid label at address ", current_address)
                    exit()
                    #return
            return imm[0] + imm[2:8] + str(rs2) + str(rs1) + str(funct3) + imm[8:12] + imm[1] + str(opcode)

        if instruction[0] in ["lw"]:
            if(len(arguments)!=3):
                print("Invalid instruction at address ", current_address)
                exit()
            
            if(arguments[0] not in registers or arguments[2] not in registers):
                print("Invalid register at address ", current_address)
                exit()
            
            rd = registers[arguments[0]]
            rs1 = registers[arguments[2]]
            
            if(int(arguments[1]) >= 2**11 or int(arguments[1]) < -2**11):
                print("Invalid immediate value")
                exit()
                #return
            
            if(int(arguments[1]) >= 0):
                imm = '{:012b}'.format(int(arguments[1]))
            else:
                imm = '{:012b}'.format(2**12 + int(arguments[1]))
            return imm + str(rs1) + str(funct3) + str(rd) + str(opcode)

        if instruction[0] in ["sw"]:
            if(len(arguments)!=3):
                print("Invalid instruction at address ", current_address)
                exit()
            
            if(arguments[0] not in registers or arguments[2] not in registers):
                print("Invalid register at address ", current_address)
                exit()
            
            rs2 = registers[arguments[0]]
            rs1 = registers[arguments[2]]
            if(int(arguments[1]) >= 2**11 or int(arguments[1]) < -2**11):
                print("Invalid immediate value")
                exit()
                #return
            
            if(int(arguments[1]) >= 0):
                imm = '{:012b}'.format(int(arguments[1]))
            else:
                imm = '{:012b}'.format(2**12 + int(arguments[1]))
            return imm[0:7] + str(rs2) + str(rs1) + str(funct3) + imm[7:] + str(opcode)

        if instruction[0] in ["addi", "sltiu"]:
            if(len(arguments)!=3):
                print("Invalid instruction at address ", current_address)
                exit()
            
            if(arguments[0] not in registers or arguments[1] not in registers):
                print("Invalid register at address ", current_address)
                exit()
            
            rd = registers[arguments[0]]
            rs1 = registers[arguments[1]]
            if(int(arguments[2]) >= 2**11 or int(arguments[2]) < -2**11):
                print("Invalid immediate value")
                exit()
                #return
            if(int(arguments[2]) >= 0):
                imm = '{:012b}'.format(int(arguments[2]))
            else:
                imm = '{:012b}'.format(2**12 + int(arguments[2]))
            if funct7 is None:
                return imm + str(rs1) + str(funct3) + str(rd) + str(opcode)
            else:
                return imm + str(rs1) + str(funct3) + str(rd) + str(opcode)

        if instruction[0] in ["add", "sub", "sll", "slt", "sltu", "xor", "srl", "or", "and"]:
            if(len(arguments)!=3):
                print("Invalid instruction at address ", current_address)
                exit()
            
            if(arguments[0] not in registers or arguments[1] not in registers or arguments[2] not in registers):
                print("Invalid register at address ", current_address)
                exit()
            
            rd = registers[arguments[0]]
            rs1 = registers[arguments[1]]
            rs2 = registers[arguments[2]]
            return str(funct7) + str(rs2) + str(rs1) + str(funct3) + str(rd) + str(opcode)
        
        
    def assemble(self):
        self.machine_code = []
        for i in range(len(self.instructions)):
            self.current_address = 4 * i
            if(self.asm1[i]!=[None]):
                assembled_instruction = self.assemble_instruction(self.asm1[i], self.current_address)
                if(assembled_instruction != None):
                    self.machine_code.append(assembled_instruction)
                else:
                    exit()
            
            print(self.current_address)

        print(self.machine_code)
        
    def write_to_file(self, filename):
        with open(filename, "w") as f:
            for i in range(len(self.machine_code)-1):
                if(self.machine_code[i] != None):
                    f.write(self.machine_code[i] + "\n")
                else:
                    f.write("\n")
            f.write(self.machine_code[-1])

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("usage: assembler.py <input_filename> <output_filename>")
        sys.exit(1)
    assembler = Assembler(sys.argv[1])
    assembler.assemble()
    assembler.write_to_file(sys.argv[2])