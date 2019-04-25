package picard.arrays.illumina;

import htsjdk.samtools.util.IOUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import org.testng.Assert;
import org.testng.annotations.Test;
import picard.PicardException;
import picard.arrays.GtcToVcf;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

public class GtcToVcfTest {
    private static final File TEST_DATA_DIR = new File("testdata/picard/arrays/illumina");
    private static final File TEST_EXTENDED_MANIFEST_FILE = new File(TEST_DATA_DIR, "Test.extended.csv");
    private static final File TEST_GTC_RECORDS_FILE = new File(TEST_DATA_DIR, "Test.gtc_records.csv");

    private static final List<String> INFINIUM_VCF_FORMAT_FIELDS = Arrays.asList(
            InfiniumVcfFields.X,
            InfiniumVcfFields.Y,
            InfiniumVcfFields.NORMX,
            InfiniumVcfFields.NORMY,
            InfiniumVcfFields.R,
            InfiniumVcfFields.THETA,
            InfiniumVcfFields.BAF,
            InfiniumVcfFields.LRR,
            InfiniumVcfFields.IGC);

    @Test
    public void testGetGenotype() throws IOException {
        GtcToVcf gtcToVcf = new GtcToVcf();
        final ExtendedIlluminaManifest manifest = new ExtendedIlluminaManifest(TEST_EXTENDED_MANIFEST_FILE);
        final List<InfiniumGTCRecord> infiniumGTCRecords = loadInfiniumGTCRecords();

        final Iterator<ExtendedIlluminaManifestRecord> iterator = manifest.extendedIterator();
        int gtcIndex = 0;
        while (iterator.hasNext()) {
            final ExtendedIlluminaManifestRecord record = iterator.next();
            final InfiniumGTCRecord infiniumGtcRecord = infiniumGTCRecords.get(gtcIndex++);
            Allele A = record.getAlleleA();
            Allele B = record.getAlleleB();

            // The Sample Alleles
            final List<Allele> alleles;

            if (infiniumGtcRecord.genotype == InfiniumGTCFile.NO_CALL) alleles = GtcToVcf.NO_CALL_ALLELES;
            else if (infiniumGtcRecord.genotype == InfiniumGTCFile.AA_CALL) alleles = Arrays.asList(A, A);
            else if (infiniumGtcRecord.genotype == InfiniumGTCFile.AB_CALL) alleles = Arrays.asList(A, B);
            else if (infiniumGtcRecord.genotype == InfiniumGTCFile.BB_CALL) alleles = Arrays.asList(B, B);
            else {
                throw new PicardException("Unexpected genotype call [" + infiniumGtcRecord.genotype + "]" + " for SNP: " + record.getName());
            }

            Genotype genotype = gtcToVcf.getGenotype("test", infiniumGtcRecord, record, A, B);
            Allele genotypeAllele1 = genotype.getAllele(0);
            Allele genotypeAllele2 = genotype.getAllele(1);
            Assert.assertTrue(genotypeAllele1.basesMatch(alleles.get(0)));
            Assert.assertTrue(genotypeAllele2.basesMatch(alleles.get(1)));

            Assert.assertEquals(genotype.getSampleName(), "test");
            Map<String, Object> extendedAttributes = genotype.getExtendedAttributes();
            Assert.assertTrue(extendedAttributes.keySet().containsAll(INFINIUM_VCF_FORMAT_FIELDS));
            Assert.assertEquals(extendedAttributes.get(InfiniumVcfFields.IGC), GtcToVcf.formatFloatForVcf(infiniumGtcRecord.genotypeScore));
            Assert.assertEquals(extendedAttributes.get(InfiniumVcfFields.X), infiniumGtcRecord.rawXIntensity);
            Assert.assertEquals(extendedAttributes.get(InfiniumVcfFields.Y), infiniumGtcRecord.rawYIntensity);
            Assert.assertEquals(extendedAttributes.get(InfiniumVcfFields.NORMX), GtcToVcf.formatFloatForVcf(infiniumGtcRecord.normalizedXIntensity));
            Assert.assertEquals(extendedAttributes.get(InfiniumVcfFields.NORMY), GtcToVcf.formatFloatForVcf(infiniumGtcRecord.normalizedYIntensity));
            Assert.assertEquals(extendedAttributes.get(InfiniumVcfFields.R), GtcToVcf.formatFloatForVcf(infiniumGtcRecord.RIlmn));
            Assert.assertEquals(extendedAttributes.get(InfiniumVcfFields.THETA), GtcToVcf.formatFloatForVcf(infiniumGtcRecord.thetaIlmn));
            Assert.assertEquals(extendedAttributes.get(InfiniumVcfFields.BAF), GtcToVcf.formatFloatForVcf(infiniumGtcRecord.bAlleleFreq));
            Assert.assertEquals(extendedAttributes.get(InfiniumVcfFields.LRR), GtcToVcf.formatFloatForVcf(infiniumGtcRecord.logRRatio));
            Assert.assertEquals(extendedAttributes.get(InfiniumVcfFields.IGC), GtcToVcf.formatFloatForVcf(infiniumGtcRecord.genotypeScore));

            // Test that there aren't any other attributes added.
            Set<String> foundKeys = extendedAttributes.keySet();
            foundKeys.removeAll(INFINIUM_VCF_FORMAT_FIELDS);
            Assert.assertTrue(foundKeys.isEmpty(), "Found extra attributes in FORMAT field of VCF");
        }
    }

    private List<InfiniumGTCRecord> loadInfiniumGTCRecords() throws FileNotFoundException {
        final List<String> lines = IOUtil.slurpLines(TEST_GTC_RECORDS_FILE);
        final List<InfiniumGTCRecord> infiniumGTCRecords = new ArrayList<>();
        for (String line : lines) {
            infiniumGTCRecords.add(new InfiniumGTCRecord(line));
        }
        return infiniumGTCRecords;
    }
}

