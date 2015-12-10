package de.dkfz.b080.co.qcworkflow;

import de.dkfz.b080.co.files.*;
import de.dkfz.roddy.execution.jobs.*;
import de.dkfz.roddy.knowledge.files.Tuple3;
import de.dkfz.roddy.knowledge.methods.GenericMethod;

/**
 * Created by kleinhei on 6/16/14.
 */
@StaticScriptProviderClass
public final class AceSeqQC {


    @ScriptCallingMethod
    public static TextFile annotateCovWindows(CoverageTextFile coverage1kbFile) {
        return (TextFile) GenericMethod.callGenericTool(COConstants.TOOL_ADD_HAPLOTYPES_TO_SNP_FILE, coverage1kbFile);
    }

    @ScriptCallingMethod
    public static TextFile mergeAndFilterCovWindows(TextFile annotatedCoverageFile) {
        return (TextFile) GenericMethod.callGenericTool(COConstants.TOOL_ADD_HAPLOTYPES_TO_SNP_FILE, annotatedCoverageFile);
    }

    @ScriptCallingMethod
    public static Tuple3<TextFile, TextFile, TextFile> correctGC(TextFile mergedAndFilteredCovWinFile) {
        return (Tuple3<TextFile, TextFile, TextFile>) GenericMethod.callGenericTool(COConstants.TOOL_CORRECT_GC_BIAS, mergedAndFilteredCovWinFile);
    }

}

