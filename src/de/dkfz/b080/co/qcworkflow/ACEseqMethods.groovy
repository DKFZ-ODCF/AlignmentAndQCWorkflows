package de.dkfz.b080.co.aceseq;

import de.dkfz.b080.co.files.*;
import de.dkfz.roddy.Roddy;
import de.dkfz.roddy.core.ExecutionContext;
import de.dkfz.roddy.core.ExecutionContextError;
import de.dkfz.roddy.execution.io.fs.FileSystemInfoProvider;
import de.dkfz.roddy.execution.jobs.*;
import de.dkfz.roddy.knowledge.files.BaseFile;
import de.dkfz.roddy.knowledge.files.Tuple2;
import de.dkfz.roddy.knowledge.files.FileGroup;
import de.dkfz.roddy.knowledge.files.Tuple3;
import de.dkfz.roddy.knowledge.methods.GenericMethod;


import javax.xml.soap.Text;
import java.io.File;
import java.nio.file.spi.FileSystemProvider;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import static de.dkfz.roddy.execution.io.fs.FileSystemInfoProvider.*;

/**
 * Created by kleinhei on 6/16/14.
 */
@StaticScriptProviderClass
public final class ACESeqMethods {


    @ScriptCallingMethod
    public static TextFile annotateCovWindows(TextFile coverage1kbFile) {
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

