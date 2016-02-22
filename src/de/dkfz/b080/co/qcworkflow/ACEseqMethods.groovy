package de.dkfz.b080.co.qcworkflow;

import de.dkfz.b080.co.files.*;
import de.dkfz.roddy.execution.jobs.*;
import de.dkfz.roddy.knowledge.files.Tuple2;
import de.dkfz.roddy.knowledge.files.Tuple3;
import de.dkfz.roddy.knowledge.methods.GenericMethod

/**
 * Created by kleinhei on 6/16/14.
 */
@StaticScriptProviderClass
public final class ACEseqMethods {


    @ScriptCallingMethod
    public static Tuple2<TextFile, TextFile> annotateCovWindows(CoverageTextFile coverage1kbFile) {
        return (Tuple2<TextFile, TextFile>) GenericMethod.callGenericTool("annotateCovWindows", coverage1kbFile);
    }

    @ScriptCallingMethod
    public static TextFile mergeAndFilterCovWindows(TextFile annotatedCoverageFile) {
        return (TextFile) GenericMethod.callGenericTool("mergeAndFilterCovWindows", annotatedCoverageFile);
    }

    @ScriptCallingMethod
    public static Tuple3<TextFile, TextFile, TextFile> correctGc(TextFile mergedAndFilteredCovWinFile) {
        return (Tuple3<TextFile, TextFile, TextFile>) GenericMethod.callGenericTool("correctGc", mergedAndFilteredCovWinFile);
    }

    private static Tuple3<TextFile, TextFile, TextFile>  aceSeqQc(CoverageTextFile windowedCoverageTextFile) {
        Tuple2<TextFile, TextFile> annotationResult = ACEseqMethods.annotateCovWindows(windowedCoverageTextFile);
        TextFile mergedAndFilteredCoverageWindowFiles = ACEseqMethods.mergeAndFilterCovWindows(annotationResult.value0);
        Tuple3<TextFile, TextFile, TextFile> correctedWindowFile = ACEseqMethods.correctGc(mergedAndFilteredCoverageWindowFiles);
        return new Tuple3(annotationResult, mergedAndFilteredCoverageWindowFiles, correctedWindowFile)
    }



}

