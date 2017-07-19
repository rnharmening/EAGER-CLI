/*
 * Copyright (c) 2016. EAGER-CLI Alexander Peltzer
 * This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package Modules.preprocessing;
import IO.Communicator;
import Modules.AModule;

import java.util.Map;

/**
 * Created by peltzer on 22.01.14.
 */
public class EAGERVersions extends AModule {

    public EAGERVersions(Communicator c) {
        super(c);
        setParameters();
    }



    @Override
    public void setProcessEnvironment (Map <String, String> env) {
        if ( !this.communicator.isUsesystemtmpdir() ) {
          AModule.setEnvironmentForParameterPrepend (env,
                                                     " ",
                                                     "JAVA_TOOL_OPTIONS",
                                                     "-Djava.io.tmpdir=" + getOutputfolder() + System.getProperty ("file.separator") + ".tmp");
        }
    }

    @Override
    public void setParameters() {

        this.parameters =  new String[]{"eagerVersions", this.communicator.getGUI_resultspath()+"/EAGER.versions"};
        this.outputfile = this.inputfile;

    }

    @Override
    public String getOutputfolder() {
        return this.communicator.getGUI_resultspath()+ "/0-FastQC";
    }


    @Override
    public String getModulename(){
        return super.getModulename();
    }



}
