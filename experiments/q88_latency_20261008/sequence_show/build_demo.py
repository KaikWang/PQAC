"""Compile/link the standalone demo using an existing CMake target's flags."""
import shlex
import subprocess
import sys
from pathlib import Path
build=Path(sys.argv[1]).resolve()
source=Path(__file__).with_name('sequence_show_demo.c').resolve()
obj=source.with_suffix('.o')
binary=source.with_suffix('')
target=build/'CMakeFiles/revocable_ac_demo.dir'
flags={}
for line in (target/'flags.make').read_text().splitlines():
 if line.startswith(('C_FLAGS =','C_INCLUDES =','C_DEFINES =')):
  key,value=line.split('=',1);flags[key.strip()]=shlex.split(value)
subprocess.run(['/usr/bin/cc',*flags['C_FLAGS'],*flags['C_DEFINES'],*flags['C_INCLUDES'],'-c',str(source),'-o',str(obj)],check=True)
link=shlex.split((target/'link.txt').read_text())
for i,value in enumerate(link):
 if value.endswith('/revocable_ac_demo.c.o'):link[i]=str(obj)
 if i>0 and link[i-1]=='-o':link[i]=str(binary)
subprocess.run(link,cwd=build,check=True)
