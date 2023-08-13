#include <aws/core/Aws.h>
#include <aws/s3/S3Client.h>
#include <aws/s3/model/CreateBucketRequest.h>

// From https://stackoverflow.com/questions/74869923/how-do-i-stream-data-from-retrieved-getobject-file-to-localfile-on-disk-c
#ifdef EXAMPKLE!
bool GetObject(const Aws::String& objectKey,
    const Aws::String& fromBucket,
    const Aws::Client::ClientConfiguration& clientConfig) {
    Aws::S3::S3Client client(clientConfig);

    Aws::S3::Model::GetObjectRequest request;
    request.SetBucket(fromBucket);
    request.SetKey(objectKey);

    Aws::S3::Model::GetObjectOutcome outcome =
        client.GetObject(request);

    if (!outcome.IsSuccess()) {
        const Aws::S3::S3Error& err = outcome.GetError();
        std::cerr << "Error: GetObject: " <<
            err.GetExceptionName() << ": " << err.GetMessage() << std::endl;
    }
    else {
        std::cout << "Successfully retrieved '" << objectKey << "' from '"
            << fromBucket << "'." << std::endl;

        //create the filename, which will be the objectKey 
        std::string local_file_name = "./netcdf/" + objectKey;
        std::ofstream local_file(local_file_name, std::ios::binary);
        auto &retrieved = outcome.GetResult().GetBody();
        local_file << retrieved.rdbuf();
        std::cout << "Done!";

    }

    return outcome.IsSuccess();
}
#endif  // EXAMPKLE


namespace aws_read_file {
}  // namespace read_file

int
main() {
}

