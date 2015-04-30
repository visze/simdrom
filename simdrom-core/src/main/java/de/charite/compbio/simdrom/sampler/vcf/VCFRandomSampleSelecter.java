/**
 * 
 */
package de.charite.compbio.simdrom.sampler.vcf;

import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;

import java.io.File;
import java.util.List;
import java.util.Random;

/**
 * @author Max Schubach <max.schubach@charite.de>
 *
 */
public class VCFRandomSampleSelecter {

	private VCFFileReader parser;
	private String sample;

	public VCFRandomSampleSelecter(String filePath) {
		this.parser = new VCFFileReader(new File(filePath), false);
	}

	public String getSample() {
		if (sample == null)
			sample = selectSample();
		return sample;
	}
	
	private String selectSample() {
		VCFHeader test = parser.getFileHeader();
		System.out.println(test);
		List<String> samples = test.getGenotypeSamples();
		parser.close();
		int num = new Random().nextInt(samples.size());
		return samples.get(num);
	}


}
